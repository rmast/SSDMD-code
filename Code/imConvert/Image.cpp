#include "include/Image.hpp"
#include <sys/stat.h>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <set>
#include "include/ImageWriter.hpp"
#include "include/skeleton_cuda.hpp"
#include "include/connected.hpp"
#include "include/messages.h"
#include <parallel/algorithm>
#include <unordered_set>
#include <boost/functional/hash.hpp>

string OUTPUT_FILE;
bool BUNDLE, OVERLAP_PRUNE, MSIL_SELECTION;
float ALPHA;
float OVERLAP_PRUNE_EPSILON;
int EPSILON;
int MIN_PATH_LENGTH;
int MIN_SUM_RADIUS;
int MIN_OBJECT_SIZE;
int peaks;
int distinguishable_interval;
#define EPSILON 0.00001

/*************** CONSTRUCTORS ***************/
Image::Image(FIELD<float> *in, float islandThresh, float importanceThresh, int GrayInterval) {
    PRINT(MSG_NORMAL, "Creating Image Object...\n");
    this->layerThreshold = importanceThresh;
    this->DistinguishableInterval = GrayInterval;
    this->islandThreshold = islandThresh;
    this->importance = NULL;
    this->im = in;
    this->nPix = in->dimX() * in->dimY();
    std::set<int> levels(in->data(), in->data() + this->nPix);
    this->graylevels = reinterpret_cast<int*>(malloc(levels.size() * sizeof(int)));
    std::set<int>::iterator it = levels.begin();
    for (uint i = 0; i < levels.size(); ++i) {
        this->graylevels[i] = *it;
        std::advance(it, 1);
    }
    this->numLayers = levels.size();
    PRINT(MSG_NORMAL, "Done!\n");
}

Image::Image(FIELD<float> *in) {
    Image(in, 0, 0, 0);
}

Image::~Image() {
    deallocateCudaMem();
    free(importance);
    free(graylevels);
}

/*************** FUNCTIONS **************/
void detect_peak(
    const double*   data, /* the data */
    int             data_count, /* row count of data */
    vector<int>&    emi_peaks, /* emission peaks will be put here */
    double          delta, /* delta used for distinguishing peaks */
    int             emi_first /* should we search emission peak first of
                                     absorption peak first? */
) {
    int     i;
    double  mx;
    double  mn;
    int     mx_pos = 0;
    int     mn_pos = 0;
    int     is_detecting_emi = emi_first;


    mx = data[0];
    mn = data[0];

    for (i = 1; i < data_count; ++i) {
        if (data[i] > mx) {
            mx_pos = i;
            mx = data[i];
        }
        if (data[i] < mn) {
            mn_pos = i;
            mn = data[i];
        }

        if (is_detecting_emi &&
                data[i] < mx - delta) {

            emi_peaks.push_back(mx_pos);
            is_detecting_emi = 0;

            i = mx_pos - 1;

            mn = data[mx_pos];
            mn_pos = mx_pos;
        } else if ((!is_detecting_emi) &&
                   data[i] > mn + delta) {

            is_detecting_emi = 1;

            i = mn_pos - 1;

            mx = data[mn_pos];
            mx_pos = mn_pos;
        }
    }
}

void find_peaks(double* importance, double width) {
    double impfac = 0.1;
    vector<int> v;
    int numiters = 0;
    while (numiters < 1000) {
        v.clear();
        detect_peak(importance, 256, v, impfac, 0);
        if (v.size() < width)
            impfac *= .9;
        else if (v.size() > width)
            impfac /= .9;
        else
            break;
        numiters++;
    }
    memset(importance, 0, 256 * sizeof(double));
    for (auto elem : v)
        importance[elem] = 1;
}

void detect_layers(double* upper_level, double threshold, bool needAssign)//float->double, is there a problem?
{
    peaks = 0;
    int i = 0;
    int StartPoint = distinguishable_interval; //There is no need to check i and i+1; it is not distinguished by eyes, so check between i and i+StartPoint.
    double difference;

    double copy_upper_level[256];
    for (int j = 0; j < 256; ++j){
        copy_upper_level [j] = upper_level [j];
    }

    while ((i + StartPoint) < 256)
    {
        difference = copy_upper_level[i] - copy_upper_level[i + StartPoint];//attention: here must not be upper_level
        if (difference > threshold)//choose this layer
        {
            if (needAssign) {
                upper_level[i + StartPoint] = 1;
            }
            
            i = i + StartPoint;
           // PRINT(MSG_NORMAL, "Choose_layers: %d\n",i);
            StartPoint = distinguishable_interval;
            peaks++;
        }
        else
        {
            StartPoint += 1;
        }
        
    }

}

//binary search
void find_layers(double* importance_upper, double width)
{
    double impfac = 0.5;
    double head = 0;
    double tail = 0.5;
    int numiters = 0;
    while (numiters < 50) {
       
        detect_layers(importance_upper, impfac, 0);
        //PRINT(MSG_NORMAL, "the number of peaks: %d\n",peaks);
        if (peaks < width){//impfac need a smaller one
            tail = impfac;
            impfac = (head + impfac)/2;
        }
            
        else if (peaks > width) //impfac need a bigger one
            {
                head = impfac;
                impfac = (tail + impfac)/2;
            }
            else
             break;
        numiters++;
    }
   // PRINT(MSG_NORMAL, "numiter:%d...\n", numiters);
  // PRINT(MSG_NORMAL, "numiter:%f...\n", impfac);
   
    detect_layers(importance_upper, impfac, 1);
    
     for (int i = 0; i < 256; ++i) {
        if (importance_upper[i] != 1) {
            importance_upper[i] = 0;
        }
    }
}

/*
void find_layers(double* importance_upper, double width)
{
    double impfac = 0.02;
    
    int numiters = 0;
    while (numiters < 1000) {
       
        detect_layers(importance_upper, impfac, 0);
        //PRINT(MSG_NORMAL, "the number of peaks: %d\n",peaks);
        if (peaks < width){
            impfac *= .9;
        }
            
        else if (peaks > width)
            impfac /= .9;
        else
            break;
        numiters++;
    }
   PRINT(MSG_NORMAL, "numiter:%d...\n", numiters);
   
    detect_layers(importance_upper, impfac, 1);
    
     for (int i = 0; i < 256; ++i) {
        if (importance_upper[i] != 1) {
            importance_upper[i] = 0;
        }
    }
}


*/

/*
* Calculate the histogram of the image, which is equal to the importance for each level.
* Avoid the use of in->value(), because it is less efficient (performs multiplications).
* The order is irrelevant anyway.
*/
void Image::calculateImportance() {
    PRINT(MSG_NORMAL, "Calculating the importance for each layer...\n");
    int normFac = 0;
    float *c = im->data();
    float *end = im->data() + nPix;
    /* If importance was already calculated before, cleanup! */
    if (importance) { free(importance); }
    importance = static_cast<double*>(calloc(256, sizeof(double)));
    double* UpperLevelSet = static_cast<double*>(calloc(256, sizeof(double)));

    if (!importance) {
        PRINT(MSG_ERROR, "Error: could not allocate importance histogram\n");
        exit(-1);
    }
    // Create a histogram
    while (c < end) {
        importance[static_cast<int>(*c++)] += 1;
    }

    ////////upper_level set histogram/////////
    UpperLevelSet[255] = importance[255];
    for (int i = 255; i > 0; i--)
    {
        UpperLevelSet[i-1] = importance[i-1] + UpperLevelSet[i];
    }

    // Normalize it
    //normFac = static_cast<int>(*std::max_element(importance, importance + 256));//find the max one.
    double max = UpperLevelSet[0];
    for (int i = 0; i < 256; ++i) {
        UpperLevelSet[i] = UpperLevelSet[i] / static_cast<double>(max) - EPSILON;//To avoid to be 1.
    }

    // Find the relevant layers based on peaks in the histogram
    if (MSIL_SELECTION) {
    
        distinguishable_interval = DistinguishableInterval;
        find_layers(UpperLevelSet, layerThreshold);
        //memcpy(importance, UpperLevelSet, sizeof(UpperLevelSet));
        importance = UpperLevelSet;
    } 
    else {
        // else local-maxima method
        find_peaks(importance, layerThreshold);
    }


    std::vector<int> v;
    for (int i = 0; i < 256; ++i) {
        if (importance[i] == 1) {
            v.push_back(i);
        }
    }
    if (v.size() == 0) {
        PRINT(MSG_ERROR, "ERROR: No layers selected. Exiting...\n");
        exit(-1);
    }
    PRINT(MSG_NORMAL, "Selected %lu layers: ", v.size());
    std::ostringstream ss;

    std::copy(v.begin(), v.end() - 1, std::ostream_iterator<int>(ss, ", "));
    ss << v.back();
    PRINT(MSG_NORMAL, "%s\n", ss.str().c_str());
    v.clear();
}

/**
* fullDilate and fullErode are placeholders. Although they do work, they
* should be replaced by better erode and dilation functions. This is used primarily
* to test how much an opening on the skeleton will reduce the image size.
*/

/* fullDilate -- Perform dilation with a S.E. of 3x3, filled with ones. */
FIELD<float> * fullDilate(FIELD<float> *layer) {
    FIELD<float> *ret = new FIELD<float>(layer->dimX(), layer->dimY());
    memset(ret->data(), 0, layer->dimX() * layer->dimY() * sizeof(float));
    for (int y = 0; y < layer->dimY(); ++y) {
        for (int x = 0; x < layer->dimX(); ++x) {
            if (layer->value(x, y)) {
                ret->set(x - 1, y - 1, 255);
                ret->set(x - 1, y  , 255);
                ret->set(x - 1, y + 1, 255);
                ret->set(x  , y - 1, 255);
                ret->set(x  , y  , 255);
                ret->set(x  , y + 1, 255);
                ret->set(x + 1, y - 1, 255);
                ret->set(x + 1, y  , 255);
                ret->set(x + 1, y + 1, 255);
            }
        }
    }
    delete layer;
    return ret;
}

/* fullErode -- Perform erosion with a S.E. of 3x3, filled with ones. */
FIELD<float> * fullErode(FIELD<float> *layer) {
    FIELD<float> *ret = new FIELD<float>(layer->dimX(), layer->dimY());
    for (int y = 0; y < layer->dimY(); ++y) {
        for (int x = 0; x < layer->dimX(); ++x) {
            if (
                layer->value(x - 1, y - 1) &&
                layer->value(x - 1, y) &&
                layer->value(x - 1, y + 1) &&
                layer->value(x, y - 1) &&
                layer->value(x, y) &&
                layer->value(x, y + 1) &&
                layer->value(x + 1, y - 1) &&
                layer->value(x + 1, y) &&
                layer->value(x + 1, y + 1)
            ) {
                ret->set(x, y, 255);
            } else { ret->set(x, y, 0); }
        }
    }
    delete layer;
    return ret;
}

/* rmObject -Remove current object in a 3x3 kernel, used for removeDoubleSkel: */
void rmObject(int *k, int x, int y) {
    if (x < 0 || x > 2 || y < 0 || y > 2 || k[y * 3 + x] == 0) { return; }
    k[y * 3 + x] = 0;
    rmObject(k, x + 1, y + 1);
    rmObject(k, x + 1, y);
    rmObject(k, x + 1, y - 1);
    rmObject(k, x, y + 1);
    rmObject(k, x, y - 1);
    rmObject(k, x - 1, y + 1);
    rmObject(k, x - 1, y);
    rmObject(k, x - 1, y - 1);
}
/* numObjects - Count the number of objects in a 3x3 kernel, used for removeDoubleSkel: */
int numObjects(int *k) {
    int c = 0;
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; ++y) {
            if (k[y * 3 + x]) { rmObject(k, x, y); c++; }
        }
    }
    return c;
}
/* End count code */

/**
* removeDoubleSkel
* @param FIELD<float> * layer -- the layer of which the skeleton should be reduced
* @return new FIELD<float> *. Copy of 'layer', where all redundant skeleton-points are removed (i.e. rows of width 2.)
*/
FIELD<float> * removeDoubleSkel(FIELD<float> *layer) {
    int *k = (int *)calloc(9, sizeof(int));
    for (int y = 0; y < layer->dimY(); ++y) {
        for (int x = 0; x < layer->dimX(); ++x) {
            if (layer->value(x, y)) {
                k[0] = layer->value(x - 1, y - 1);
                k[1] = layer->value(x - 1, y);
                k[2] = layer->value(x - 1, y + 1);
                k[3] = layer->value(x  , y - 1);
                k[4] = 0;
                k[5] = layer->value(x  , y + 1);
                k[6] = layer->value(x + 1, y - 1);
                k[7] = layer->value(x + 1, y);
                k[8] = layer->value(x + 1, y + 1);
                if (k[0] + k[1] + k[2] + k[3] + k[4] + k[5] + k[6] + k[7] + k[8] > 256) {
                    int b = numObjects(k);
                    if (b < 2 ) {layer->set(x, y, 0); }
                }
            }
        }
    }
    free(k);
    return layer;
}


/**
* Given a binary layer, remove all islands smaller than iThresh.
* , where k is the current intensity.
* @param layer
* @param iThresh
*/
void Image::removeIslands(FIELD<float>*layer, unsigned int iThresh) {
    int                     nPix    = layer->dimX() * layer->dimY();
    ConnectedComponents     *CC     = new ConnectedComponents(255);
    int                     *labeling = new int[nPix];
    float                   *fdata  = layer->data();
    int                     highestLabel;
    unsigned int            *num_pix_per_label;

    /* CCA -- store highest label in 'max' -- Calculate histogram */
    highestLabel = CC->connected(fdata, labeling, layer->dimX(), layer->dimY(), std::equal_to<float>(), true);
    // Count the number of pixels for each label
    num_pix_per_label = (unsigned int*) calloc(highestLabel + 1, sizeof(unsigned int));
    for (int j = 0; j < nPix; j++) { num_pix_per_label[labeling[j]]++; }

    /* Remove small islands based on a percentage of the image size*/
    for (int j = 0; j < nPix; j++) {
        fdata[j] = (num_pix_per_label[labeling[j]] >= (iThresh / 100.0) * nPix) ? fdata[j] : 255 - fdata[j];
    }

    /* Cleanup */
    free(num_pix_per_label);
    delete [] labeling;
    delete CC;
}


/*
* Remove small islands according the the islandThreshold variable. Notice that both "on" and "off"
* regions will be filtered.
*/
void Image::removeIslands(FIELD<float>* sm, float flagThreshold ) {
    int i, j, k;                    /* Iterators */
    FIELD<float> *inDuplicate = 0;  /* Duplicate, because of inplace modifications */
    FIELD<float> *newImg = new FIELD<float>(im->dimX(), im->dimY());
    FIELD<float> *LabelImg = new FIELD<float>(im->dimX(), im->dimY());
    int highestLabel;               /* for the CCA */
    int *ccaOut;                    /* labeled output */
    ConnectedComponents *CC;        /* CCA-object */
    float *fdata;
    unsigned int *hist;
    float *flag;//wang
    float* smd = sm -> data();//saliency map data. wang
    //for (int j = 0; j < nPix; j+=10) {
   // PRINT(MSG_NORMAL, "smd[j]: %f\n",smd[j]);} //check data

    PRINT(MSG_NORMAL, "Removing small islands...\n");
    /* Some basic initialization */
    memset(newImg->data(), 0, nPix * sizeof(float));
    memset(LabelImg->data(), 0, nPix * sizeof(float));

    /* Connected Component Analysis */
    #pragma omp parallel for private(i, j, k, ccaOut, CC, fdata, highestLabel, hist, inDuplicate)

    for (i = 0; i < 0xff; ++i) {
        PRINT(MSG_VERBOSE, "Layer: %d\n", i);
        // The below value refers to the expected number of components in an image.
        CC = new ConnectedComponents(255);
        ccaOut = new int[nPix];

        inDuplicate = (*im).dupe();
        inDuplicate->threshold(i);//threshold-set..
        
       // if (i==95)
        //   inDuplicate->writePGM("gray95.pgm");//debug
       
       
        fdata = inDuplicate->data();

        /* CCA -- store highest label in 'max' -- Calculate histogram */
        highestLabel = CC->connected(fdata, ccaOut, im->dimX(), im->dimY(), std::equal_to<float>(), true);
        hist = static_cast<unsigned int*>(calloc(highestLabel + 1, sizeof(unsigned int)));
        if (!hist) {
            PRINT(MSG_ERROR, "Error: Could not allocate histogram for connected components\n");
            exit(-1);
        } 
        //wang
        flag = static_cast<float*>(calloc(highestLabel + 1, sizeof(float)));
        if (!flag) {
            PRINT(MSG_ERROR, "Error: Could not allocate for flag\n");
            exit(-1);
        }
        for (int f = 0; f < highestLabel+1; f++) { 
            flag[f] = 0.0; //init
            }
        //for (j = 0; j < nPix; j++) { hist[ccaOut[j]]++; }
        //wang
        for (j = 0; j < nPix; j++) {
            hist[ccaOut[j]] ++; 
            flag[ccaOut[j]] += (smd[j]/255.0);
        }
       // if (i==100)
       // for (int f = 0; f < highestLabel+1; f++) { 
       // PRINT(MSG_NORMAL, "flag[%d]: %f\n",f,flag[f]);}

        /* Remove small islands */
       if (flagThreshold == 0.0){//use original method
           for (j = 0; j < nPix; j++) 
            fdata[j] = (hist[ccaOut[j]] >= (islandThreshold/100*im->dimX()*im->dimY())) ? fdata[j] : 255 - fdata[j]; //change the absolute num. to %   
       }
       else//use saliency-based island removal method
       {
           for (j = 0; j < nPix; j++) 
            fdata[j] = (flag[ccaOut[j]] >= flagThreshold) ? fdata[j] : 255 - fdata[j]; //change the absolute num. to %   
       }
       
       
      // if (i==100)
       //    inDuplicate->writePGM("gray100after.pgm"); //debug
       
        
        #pragma omp critical
        {
            for (j = 0; j < im->dimY(); j++)
                for (k = 0; k < im->dimX(); k++)
                    if (0 == fdata[j * im->dimX() + k] && newImg->value(k, j) < i) { newImg->set(k, j, i); }
        }
/*
 //wang
    if (i == 200)
    {
        float Max = 0.0;
        int *newCCA = new int[nPix];
        for (j = 0; j < nPix; j++) Max = (ccaOut[j]>Max) ? ccaOut[j] : Max;//find Max;
        for (j = 0; j < nPix; j++) newCCA[j] = (ccaOut[j]/Max*255 > 255) ? 255 : ((ccaOut[j]/Max*255 < 0) ? 0: ccaOut[j]/Max*255);

        for (j = 0; j < im->dimY(); j++)
         for (k = 0; k < im->dimX(); k++)
            { LabelImg->set(k, j, newCCA[j * im->dimX() + k]); }
        LabelImg -> writePGM("LabelImg.pgm");
    }
//end
*/
        /* Cleanup */
        free(hist);
        delete [] ccaOut;
        delete CC;
        delete inDuplicate;
    }
    

    for (j = 0; j < im->dimY(); j++)
        for (k = 0; k < im->dimX(); k++)
            im->set(k, j, newImg->value(k, j));
    im -> writePGM("afterIsland.pgm");
    delete newImg;
    PRINT(MSG_NORMAL, "Done!\n");
}

/**
* Remove unimportant layers -- Filter all layers for which their importance is lower than layerThreshold
*/
void Image::removeLayers() {
    float val_up, val_dn;

    PRINT(MSG_NORMAL, "Filtering image layers...\n");
    PRINT(MSG_VERBOSE, "The following grayscale intensities are removed:\n");
    if (MSG_LEVEL == MSG_VERBOSE)
        for (int i = 0; i < 256; ++i) {
            if (importance[i] < 1) {
                PRINT(MSG_VERBOSE, "(%d, %6.5f)\n", i, importance[i]);
            }
        }


    for (int y = 0; y < im->dimY(); y++) {
        for (int x = 0; x < im->dimX(); x++) {
            val_up = im->value(x, y);
            val_dn = im->value(x, y);
            if (importance[(int)im->value(x, y)] == 1)
                continue;
            while (val_dn >= 0 || val_up <= 256) {
                if (val_dn >= 0) {
                    if (importance[(int)val_dn] == 1) {
                        im->set(x, y, val_dn);
                        break;
                    }
                }
                if (val_up <= 256) {
                    if (importance[(int)val_up] == 1) {
                        im->set(x, y, val_up);
                        break;
                    }
                }
                val_dn--;
                val_up++;
                // if (val_up < 256) {
                // }
            }
        }
    }
}

void draw_path(int intensity, FIELD<float>* image, skel_tree_t *st, uint16_t pLength, bool rightMost) {
    /* Not a leaf, does it have exactly 1 child (is it a continuous path?) ?*/
    if (st->numChildren() == 0) {
        Triple<int, int, int> p = st->getValue();
        image->set(p.first, p.second, intensity);
    }

    if (st->numChildren() == 1) {
        Triple<int, int, int> p = st->getValue();
        image->set(p.first, p.second, intensity);
        // add_to_hist(st->getValue(), st->getChild(0)->getValue());
        draw_path(intensity, image, st->getChild(0), pLength + 1, rightMost);
        return;
    }

    /* Fork coming up! */
    if (st->numChildren() > 1) {
        /* All "non-rightmost" children: */
        auto cur = st->getValue();
        image->set(cur.first, cur.second, intensity);

        for (int i = 0; i < (st->numChildren() - 1) ; ++i) {
            draw_path(intensity, image, st->getChild(i), 1, false);
        }
        /* Treat rightmost child different, pass a longer path length, so it jumps back further after being done with the last branch. */
        // add_to_hist(cur, st->getChild(st->numChildren() - 1)->getValue());
        draw_path(intensity, image, st->getChild(st->numChildren() - 1), 1 + pLength, rightMost);
    }
}


void draw_layer(int intensity, skel_tree_t* trees, FIELD<float>* image) {
    /* Remove overhead for empty layers. */
    if (trees->numChildren() == 0) { return; }

    /* Top layer are disjunct paths. Treat them as separate objects */
    for (int child = 0; child < trees->numChildren(); ++child) {
        skel_tree_t *curnode = (*trees)[child];
        draw_path(intensity, image, curnode, /*pLength = */1, /*rightMost = */true);
    }
}

void draw_skeletons(int szx, int szy, vector<std::pair<int, skel_tree_t*>>* forest) {
    FIELD<float>   *overlapped = new FIELD<float>(szx, szy);
    for (auto elem = forest->begin(); elem != forest->end(); ++elem) {
        int intensity = elem->first;
        skel_tree_t* trees = elem->second;
        draw_layer(intensity, trees, overlapped);
    }
    overlapped->writePGM("skeletons.pgm");
}

coord2D_list_t * neighbours(int x, int y, FIELD<float> *skel) {
    coord2D_list_t *neigh = new coord2D_list_t();
    int n[8] = {1, 1, 1, 1, 1, 1, 1, 1};

    /* Check if we are hitting a boundary on the image */
    if (x <= 0 )             {        n[0] = 0;        n[3] = 0;        n[5] = 0;    }
    if (x >= skel->dimX() - 1) {        n[2] = 0;        n[4] = 0;        n[7] = 0;    }
    if (y <= 0)              {        n[0] = 0;        n[1] = 0;        n[2] = 0;    }
    if (y >= skel->dimY() - 1) {        n[5] = 0;        n[6] = 0;        n[7] = 0;    }

    /* For all valid coordinates in the 3x3 region: check for neighbours*/
    if ((n[0] != 0) && (skel->value(x - 1, y - 1) > 0)) { neigh->push_back(coord2D_t(x - 1, y - 1)); }
    if ((n[1] != 0) && (skel->value(x    , y - 1) > 0)) { neigh->push_back(coord2D_t(x    , y - 1)); }
    if ((n[2] != 0) && (skel->value(x + 1, y - 1) > 0)) { neigh->push_back(coord2D_t(x + 1, y - 1)); }
    if ((n[3] != 0) && (skel->value(x - 1, y    ) > 0)) { neigh->push_back(coord2D_t(x - 1, y    )); }
    if ((n[4] != 0) && (skel->value(x + 1, y    ) > 0)) { neigh->push_back(coord2D_t(x + 1, y    )); }
    if ((n[5] != 0) && (skel->value(x - 1, y + 1) > 0)) { neigh->push_back(coord2D_t(x - 1, y + 1)); }
    if ((n[6] != 0) && (skel->value(x    , y + 1) > 0)) { neigh->push_back(coord2D_t(x    , y + 1)); }
    if ((n[7] != 0) && (skel->value(x + 1, y + 1) > 0)) { neigh->push_back(coord2D_t(x + 1, y + 1)); }

    return neigh;
}

skel_tree_t *tracePath(int x, int y, FIELD<float> *skel, FIELD<float> *dt) {
    coord2D_t n;
    coord2D_list_t *neigh;
    skel_tree_t *path;
    if (skel->value(x, y) == 0) {
        PRINT(MSG_ERROR, "Reached invalid point.\n");
        return NULL;
    }

    /* Create new node and add to root */
    path = new Node<coord3D_t>(coord3D_t(x, y, dt->value(x, y)));
    skel->set(x, y, 0);

    neigh = neighbours(x, y, skel);
    /* Add children */
    while (neigh->size() > 0) {
        n = *(neigh->begin());
        path->addChild(tracePath(n.first, n.second, skel, dt));
        delete neigh;
        neigh = neighbours(x, y, skel);
    }

    delete neigh;
    return path;
}

skel_tree_t* traceLayer(FIELD<float> *skel, FIELD<float> *dt) {
    skel_tree_t *root;
    coord3D_t rootCoord = coord3D_t(-1, -1, -1);
    root = new skel_tree_t( rootCoord );
    for (int y = 0; y < skel->dimY(); ++y) {
        for (int x = 0; x < skel->dimX(); ++x) {
            if (skel->value(x, y) > 0) {
                root->addChild(tracePath(x, y, skel, dt));
            }
        }
    }
    return root;
}

void bundle(FIELD<float>* skelCurr, FIELD<float>* currDT, FIELD<float>* prevDT, short* prev_skel_ft, int fboSize) {
    ALPHA = max(0, min(1, ALPHA));
    for (int i = 0; i < skelCurr->dimX(); ++i) {
        for (int j = 0; j < skelCurr->dimY(); ++j) {
            if (skelCurr->value(i, j) > 0) {
                // We might bundle this point
                int id = j * fboSize + i;

                // Find closest point
                int closest_x = prev_skel_ft[2 * id];
                int closest_y = prev_skel_ft[2 * id + 1];
                if (closest_x == -1 && closest_y == 0)
                    continue;
                // If it is already on a previous skeleton point, were done
                double distance = sqrt((i - closest_x) * (i - closest_x) + (j - closest_y) * (j - closest_y));
                if (distance == 0.0)
                    continue;

                // Create vector from the current point to the closest prev skel point
                double xvec = (closest_x - i);
                double yvec = (closest_y - j);

                // What would be the location with no bound on the shift?
                int trans_vec_x = ALPHA * xvec;
                int trans_vec_y = ALPHA * yvec;
                double trans_vec_len = sqrt(trans_vec_x * trans_vec_x + trans_vec_y * trans_vec_y);

                // If we don't move it at all were done
                if (trans_vec_len == 0.0)
                    continue;

                // Find the limited length of this vector
                double desired_length = min(trans_vec_len, EPSILON);
                // Rescale the vector such that it is bounded by epsilon
                int new_x_coor = round(i + (trans_vec_x / trans_vec_len) * desired_length);
                int new_y_coor = round(j + (trans_vec_y / trans_vec_len) * desired_length);
                // Verify that the new distance is at most near EPSILON (i.e. FLOOR(new_distance) <= EPSILON)
                double new_distance = sqrt((i - new_x_coor) * (i - new_x_coor) + (j - new_y_coor) * (j - new_y_coor));
                PRINT(MSG_VERBOSE, "%d: (%d, %d) -> (%d, %d) (%lf) -> (%d, %d) (%lf)\n", EPSILON, i, j, closest_x, closest_y, distance, new_x_coor, new_y_coor, new_distance);

                // compute a new radius
                double new_r = min(currDT->value(i, j), prevDT->value(new_x_coor, new_y_coor));
                // double new_r = currDT->value(i, j);
                skelCurr->set(i, j, 0);
                skelCurr->set(new_x_coor, new_y_coor, 255);
                currDT->set(new_x_coor, new_y_coor, new_r);
            }
        }
    }
}

double overlap_prune(FIELD<float>* skelCurr, FIELD<float>* skelPrev, FIELD<float>* currDT, FIELD<float>* prevDT) {
    // Bigger at same location
    int count = 0, total_skel_points = 0;
    for (int i = 0; i < skelPrev->dimX(); ++i) {
        for (int j = 0; j < skelPrev->dimY(); ++j) {
            // If there is a skeleton point...
            if (skelCurr->value(i, j) > 0 && skelPrev->value(i, j) > 0) {
                // cout << skelCurr->value(i, j) << ", " << skelPrev->value(i, j) << endl;
                total_skel_points++;
                // ...and the radius difference w.r.t. the next layer is small enough...
                if ((currDT->value(i, j)) >= (prevDT->value(i, j))) {
                    // ...we can delete the skeleton point
                    skelPrev->set(i, j, 0);
                    count++;
                }
            }
        }
    }
// Return overlap prune ratio for statistics
    if (total_skel_points == 0)
        return 0;
    double overlap_prune_ratio = static_cast<double>(count) / total_skel_points * 100.0;
    (MSG_VERBOSE, "Overlap pruned %d points out %d points (%.2f%%)\n", count, total_skel_points, overlap_prune_ratio);

    return overlap_prune_ratio;
}

FIELD<float>* perform_skeletonization(int level, FIELD<float>* upper_level_set, FIELD<float>* sm, float smThreshold) {
   
    auto skel_curr = computeSkeleton(level, upper_level_set, sm, smThreshold);

    upper_level_set->writePGM(string("temp_images/t-" + to_string(level) + ".pgm").c_str());
 
    skel_curr->writePGM(string("temp_images/skel-" + to_string(level) + ".pgm").c_str());
     
    upper_level_set->writePGM(string("temp_images/dt-" + to_string(level) + ".pgm").c_str());
  
    //auto skelft = skelft_to_field();  ////wang // To avoid the error in this function. ox = -1;//furthermore, I do not know what is the use of this line.
   
    //if (skelft)
    //    skelft->writePGM(string("temp_images/skelft-" + to_string(level) + ".pgm").c_str());
   
    return skel_curr;
}

void removeSmallPaths(skel_tree_t *st) {
    if (st->numChildren() > 1) {
        for (int i = 0; i < st->numChildren() && st->numChildren() > 1; ++i) {
            if (st->getChild(i)->numRChildren() < MIN_PATH_LENGTH || st->getChild(i)->importance() < MIN_SUM_RADIUS) {
                st->removeRChild(i);
                --i;
            }
        }
    }

    for (int i = 0; i < st->numChildren(); ++i) {
        removeSmallPaths(st->getChild(i));
    }
}


void removeSmallObjects(skel_tree_t *st) {
    for (int i = 0; i < st->numChildren(); ++i) {
        if (st->getChild(i)->numRChildren() < MIN_OBJECT_SIZE || st->getChild(i)->importance() < MIN_SUM_RADIUS) {
            st->removeRChild(i);
            --i;
        }
    }
}

/* Due to some preprocessing steps which simplify the skeleton we could generate "new points", Points which are
 actually outside of an object may now become skeleton points (e.g. by dilation). These points have a radius of
 0. If these points are at the end of a skeleton path, we can safely remove them.*/
void removeInvisibleEndpoints(skel_tree_t *st) {
    for (int i = 0; i < st->numChildren(); ++i) {
        removeInvisibleEndpoints(st->getChild(i));

        if (st->getChild(i)->numRChildren() == 0 && st->getChild(i)->getValue().third == 0) {
            st->removeRChild(i);
            --i;
        }
    }
}


void filterTree(skel_tree_t *st) {
    removeInvisibleEndpoints(st);
    // We can never delete if they're both zero, so we can safely return then
    if (!(MIN_PATH_LENGTH == 0 && MIN_SUM_RADIUS == 0))
        removeSmallPaths(st);
    // We can never delete if they're both zero, so we can safely return then
    if (!(MIN_OBJECT_SIZE == 0 && MIN_SUM_RADIUS == 0))
        removeSmallObjects(st);
}

void trace_layer(vector<std::pair<int, skel_tree_t*>>* forest, int level, FIELD<float>* skel, FIELD<float>* dt) {
    skel_tree_t* paths = traceLayer(skel, dt);
           
    forest->push_back(make_pair(level, paths));
}

void flatten_tree(skel_tree_t* tree, unordered_set<pair<int, int>, boost::hash< std::pair<int, int>>>& flat_map) {
    if (tree == nullptr || tree->numChildren() == 0)
        return;
    auto val = tree->getValue();
    flat_map.insert(make_pair(val.first, val.second));
    for (int i = 0; i < tree->numChildren(); ++i)
        flatten_tree(tree->getChild(i), flat_map);
}

void post_process(FIELD<float>* skel_curr, FIELD<float>* skel_prev, FIELD<float>* currDT, FIELD<float>* prevDT, vector<double>& avg_overlap_prune) {
    skel_curr = removeDoubleSkel(skel_curr);
    auto skel_prev_dupe = skel_prev->dupe();

    skel_tree_t* forest = traceLayer(skel_prev_dupe, prevDT);
    delete skel_prev_dupe;

    // TODO: This procedure of filtering is horribly slow, but necessary
    // Preferably this ought to be fixed in some way...
    filterTree(forest);
    unordered_set<pair<int, int>, boost::hash< pair<int, int>>> skel_pixels;

    // Delete all skeleton points from the image that are too short
    flatten_tree(forest, skel_pixels);
    // Delete the forest
    delete forest;
    for (int i = 0; i < skel_prev->dimX(); ++i) {
        for (int j = 0; j < skel_prev->dimY(); ++j) {
            if (skel_pixels.find(make_pair(i, j)) == skel_pixels.end()) {
                skel_prev->set(i, j, 0);
            }
        }
    }

    // Overlap prune if necessary
    if (OVERLAP_PRUNE) {
        avg_overlap_prune.push_back(overlap_prune(skel_curr, skel_prev, currDT, prevDT));
    }
}

void cool_effects(FIELD<float>* skel_curr, FIELD<float>* currDT, FIELD<float>* prevDT, short * prev_skel_ft, int fboSize) {
    if (BUNDLE) {
        if (prev_skel_ft) {
            bundle(skel_curr, currDT, prevDT, prev_skel_ft, fboSize);
        }
    }
}

/**
* Calculate the skeleton for each layer.
*/
vector<std::pair<int, skel_tree_t*>>* Image::computeSkeletons(int *firstLayer, FIELD<float>* sm, float smThreshold) {
    // Does nothing with CPU skeletons but sets up CUDA skeletonization
    int fboSize = initialize_skeletonization(im);
    short        *prev_skel_ft = 0;
    short        *curr_skel_ft = 0;
    FIELD<float> *imDupeCurr = 0;
    FIELD<float> *imDupePrev = 0;
    FIELD<float> *skelCurr = 0;
    FIELD<float> *skelPrev = 0;
    vector<std::pair<int, skel_tree_t*>>* forest = new vector<std::pair<int, skel_tree_t*>>();
    PRINT(MSG_NORMAL, "Computing the skeleton for all layers...\n");
    imDupePrev = im->dupe();
    imDupePrev->threshold(0);   
    skelPrev = computeSkeleton(0, imDupePrev, sm, smThreshold);
   
    prev_skel_ft = get_current_skel_ft();

    vector<double> avg_overlap_prune;
 
    int last_level = -1;
    for (int i = 1; i < 256; ++i) {
        //PRINT(MSG_NORMAL, "Layer: %d\n", i);
        if (importance[i] > 0) {
            if (last_level == -1)//first layer! Do not need to do skeletonization. But write the first layer num into Header.
                {  *firstLayer = i; 
                    last_level = 0; 
                }
            else{
            // Threshold the image
            imDupeCurr = im->dupe();
            imDupeCurr->threshold(i);
            // skeletonize into skelcurr. The upper-level-set imDupeCurr
            // magically transforms into the DT map
            
            //string s = "beforeskel" + std::to_string(i) + ".pgm";
            //imDupeCurr-> writePGM(s.c_str());

            
            skelCurr = perform_skeletonization(i, imDupeCurr, sm, smThreshold);//imDupeCurr now is turned to be DT.
          
            //PRINT(MSG_NORMAL, "Layer: %d\n", i);
            curr_skel_ft = get_current_skel_ft();

            // Perform cool effects by bundling
            skelCurr->writePGM(string(string("temp_images/skel-") + to_string(i) + string("a.pgm")).c_str());
            //if (last_level != -1)
            //    cool_effects(skelCurr, imDupeCurr, imDupePrev, prev_skel_ft, fboSize);
            // Remove extraneous points
            skelCurr->writePGM(string(string("temp_images/skel-") + to_string(i) + string("b.pgm")).c_str());
            post_process(skelCurr, skelPrev, imDupeCurr, imDupePrev, avg_overlap_prune);

        //multi-output.
            //s = "skel" + std::to_string(i) + ".pgm";
            //skelPrev-> writePGM(s.c_str());


            // Trace the tree and store it in the forest
            if (last_level != 0) 
                trace_layer(forest, last_level, skelPrev, imDupePrev);
           
            // Store the last level for the final layer
            last_level = i;

            delete skelPrev;
            delete imDupePrev;
            free(prev_skel_ft);
            // Move on to the next level by storing the current level as the
            // previous
            prev_skel_ft = curr_skel_ft;
            skelPrev = skelCurr;
            imDupePrev = imDupeCurr;
        }
        }
    }
    // Store the final layer
    trace_layer(forest, last_level, skelPrev, imDupePrev);
    delete skelCurr;
    delete imDupeCurr;
    free(curr_skel_ft);
    PRINT(MSG_NORMAL, "\nDone!\n");
    if (OVERLAP_PRUNE) {
        double sum = std::accumulate(avg_overlap_prune.begin(), avg_overlap_prune.end(), 0.0);
        double mean = sum / avg_overlap_prune.size();
        PRINT(MSG_NORMAL, "\nOverlap pruned %.2lf%% of all skeleton points.\n", mean);
    }
    return forest;
}
