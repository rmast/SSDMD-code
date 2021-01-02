/* main.cpp */

/**
 * Author   :   Yuri Meiburg / Maarten Terpstra / Jieying Wang
 * 2011 / 2016 / 2019
 * Rijksuniversiteit Groningen
 * Computational Science and Visualization

 * imConvert
 *
 * This program can read PGM/PPM files, and convert them to the SIR method,
 * which stands for "Skeletonal Image Representation". It is in a very
 * early research state, so no GUI, no broad image support, and not
 * heavily optimized for speed.
 *
 * To view SIR-files, use imShow.
 *
 */



#include <math.h>
#include "fileio/fileio.hpp"
#include <string>
#include <sys/resource.h>
#include "main.hpp"
#include <omp.h>
#include <vector>
#include "include/Image.hpp"
#include "include/io.hpp"
#include "include/messages.h"
#include "include/image.h"
#include "include/ImageWriter.hpp"
#include "configParser/include/Config.hpp"
#include <chrono>
#include <boost/algorithm/string.hpp>
using namespace std;


/* All properties that can be set in the config files: */
int MSG_LEVEL = MSG_NORMAL;
string filename_stdstr;
string filename_sm;
string compress_method;
/* Image space variables */
float islandThreshold = 0;
float layerThreshold = 0;
float flagThreshold = 0.0;//wang.
int DistinguishableInterval = 0;
float smThreshold = 1.0;
COLORSPACE c_space;
extern float SKELETON_DT_THRESHOLD;
extern float SKELETON_SALIENCY_THRESHOLD;
extern float SKELETON_ISLAND_THRESHOLD;
extern string OUTPUT_FILE;
extern string COMP_METHOD_NAME;
extern string ENCODING;
extern bool MSIL_SELECTION;
// Overlap pruning
extern bool OVERLAP_PRUNE;
extern float OVERLAP_PRUNE_EPSILON;
// Bundling variables
extern bool BUNDLE;
extern int EPSILON;
extern float ALPHA;
extern int COMP_LEVEL;
/* Tree variables */
extern int MIN_OBJECT_SIZE;
extern int MIN_SUM_RADIUS;
extern int MIN_PATH_LENGTH;
extern int time1;

FIELD<float>* SMfield;

// Set special colorspace in case of color image
// Defaults to RGB when an color image is encountered, otherwise to gray
COLORSPACE set_color_space(string c_space) {
    if (boost::iequals(c_space, "hsv")) {
        return COLORSPACE::HSV;
    } else if (boost::iequals(c_space, "ycbcr") || boost::iequals(c_space, "yuv")) {
        return COLORSPACE::YCC;
    } else if (boost::iequals(c_space, "rgb")) {
        return COLORSPACE::RGB;
    } else {
        return COLORSPACE::NONE;
    }
}
/* Set all parameters according to the configuration file parsed
 * with the Config object */
void setConfig(Config c) {
    if (!c.exists("filename")) {
        cout << "Please specify a filename" << endl;
        exit(-1);
    }
    filename_stdstr = c.pString("filename");//this is where you get the filename from the TestConfig.txt
    filename_sm = c.pString("filenameSM");//this is where you get the filename from the TestConfig.txt
    if (!c.exists("outputLevel")) {
        MSG_LEVEL = MSG_NORMAL;
    } else {
        switch (c.pString("outputLevel")[0]) {
        case 'q':
            MSG_LEVEL = MSG_NEVER;
            break;
        case 'e':
            MSG_LEVEL = MSG_ERROR;
            break;
        case 'n':
            MSG_LEVEL = MSG_NORMAL;
            break;
        case 'v':
            MSG_LEVEL = MSG_VERBOSE;
            break;
        }
    }

    // MSG_LEVEL = MSG_VERBOSE;
    DistinguishableInterval = c.exists("distinguishable_interval") ? c.pInt("distinguishable_interval") : 5;//Wang
    flagThreshold = c.exists("flagThreshold") ? c.pDouble("flagThreshold") : 1;//Wang
    smThreshold = c.exists("smThreshold") ? c.pInt("smThreshold") : 3;//Wang
    layerThreshold = c.exists("num_layers") ? c.pDouble("num_layers") : (c.exists("lThreshold") ? c.pDouble("lThreshold") : 0);
    islandThreshold = c.exists("islandThreshold") ? c.pDouble("islandThreshold") : 0;
    COMP_METHOD_NAME = c.exists("compression_method") ? c.pString("compression_method") : "lzma";
    COMP_LEVEL = c.exists("compression_level") ? c.pInt("compression") : -1;
    ENCODING = c.exists("encoding") ? c.pString("encoding") : "standard";
    //MSIL_SELECTION = c.exists("layer_selection") ? c.pString("layer_selection") != "threshold" : false;
    MSIL_SELECTION = c.exists("layer_selection") ? true : false;
    SKELETON_DT_THRESHOLD = c.exists("sdtThreshold") ? c.pDouble("sdtThreshold") : 5;
    SKELETON_SALIENCY_THRESHOLD = c.exists("ssThreshold") ? c.pDouble("ssThreshold") : 5;
    SKELETON_ISLAND_THRESHOLD = c.exists("SkeletonThreshold") ? c.pDouble("SkeletonThreshold") : 10;

    MIN_OBJECT_SIZE = c.exists("minObjSize") ? c.pInt("minObjSize") : 5;
    MIN_SUM_RADIUS = c.exists("minSumRadius") ? c.pInt("minSumRadius") : 15;
    MIN_PATH_LENGTH = c.exists("minPathLength") ? c.pInt("minPathLength") : 3;

    OVERLAP_PRUNE = c.exists("overlap_pruning") ? c.pBool("overlap_pruning") : false;
    OVERLAP_PRUNE_EPSILON = c.exists("overlap_pruning_epsilon") ? c.pDouble("overlap_pruning_epsilon") : 0.0001;
    BUNDLE = c.exists("bundle") ? c.pBool("bundle") : false;
    ALPHA = c.exists("alpha") ? min(1, max(0, c.pDouble("alpha"))) : 0;
    EPSILON = c.exists("epsilon") ? c.pInt("epsilon") : 0;
    OUTPUT_FILE = c.exists("outputFile") ? c.pString("outputFile") : "default.sir";
    c_space = c.exists("colorspace") ? set_color_space(c.pString("colorspace")) : COLORSPACE::NONE;

    /* output all settings, so we can check if it parsed correctly. 
    PRINT(MSG_NORMAL, "\n\n-------- Begin Configuration --------\n\n");
    PRINT(MSG_NORMAL, "Output level : % s\n", MSG_LEVEL == MSG_NORMAL ? "Normal" : MSG_LEVEL == MSG_VERBOSE ? "Verbose" : "None");
    PRINT(MSG_NORMAL, "File in : % s\n", filename_stdstr.c_str());
    PRINT(MSG_NORMAL, "File out : % s\n", OUTPUT_FILE.c_str());
    PRINT(MSG_NORMAL, "Small Island Removal Threshold : % f\n", islandThreshold);
    PRINT(MSG_NORMAL, "Layer Selection Method : % s\n", MSIL_SELECTION ? "MSIL" : "Threshold");
    if (MSIL_SELECTION) {
        PRINT(MSG_NORMAL, "Number of layers : % f\n", layerThreshold);
    } else {
        PRINT(MSG_NORMAL, "Layer Importance Threshold : % f\n", layerThreshold);
    }
    PRINT(MSG_NORMAL, "Skeleton DT Threshold : % f\n", SKELETON_DT_THRESHOLD);
    PRINT(MSG_NORMAL, "Skeleton Saliency Threshold : % f\n", SKELETON_SALIENCY_THRESHOLD);
    PRINT(MSG_NORMAL, "Skeleton CCA Threshold : % f\n", SKELETON_ISLAND_THRESHOLD);
    PRINT(MSG_NORMAL, "Minimum object size : % d\n", MIN_OBJECT_SIZE);
    PRINT(MSG_NORMAL, "Threshold for sum of radii : % d\n", MIN_SUM_RADIUS);
    PRINT(MSG_NORMAL, "Minimal branch length : % d\n", MIN_PATH_LENGTH);
    PRINT(MSG_NORMAL, "Compression method : % s\n", COMP_METHOD_NAME.c_str());
    PRINT(MSG_NORMAL, "Compression level : % d\n", COMP_LEVEL);
    PRINT(MSG_NORMAL, "Encoding method : % s\n", ENCODING.c_str());
    PRINT(MSG_NORMAL, "Perform Overlap Pruning : % s\n", OVERLAP_PRUNE ? "true" : "false");
    if (OVERLAP_PRUNE)
        PRINT(MSG_NORMAL, "Overlap Pruning epsilon : % f\n", OVERLAP_PRUNE_EPSILON);

    PRINT(MSG_NORMAL, "Perform bundling : % s\n", BUNDLE ? "true" : "false");
    if (BUNDLE) {
        PRINT(MSG_NORMAL, "Alpha : % f\n", ALPHA);
        PRINT(MSG_NORMAL, "Epsilon : % d\n", EPSILON);
    }
    PRINT(MSG_NORMAL, "Colorspace: %d\n", c_space);
    PRINT(MSG_NORMAL, "\n\n-------- End of Configuration --------\n\n");
    */
}

vector<std::pair<int, skel_tree_t*>>* execute_skeleton_pipeline(FIELD<float>* im, int *firstLayer) {
    // IS_init("stats.m");
    Image* il = new Image(im, islandThreshold, layerThreshold, DistinguishableInterval);
    // IS_analyseImage("raw", il->im);
    il->removeIslands(SMfield, flagThreshold);
    // IS_analyseImage("removed_islands", il->im);
    il->calculateImportance();
    il->removeLayers();
    // IS_analyseImage("removed_layers", il->im);
    auto begin = std::chrono::steady_clock::now();
    vector<std::pair<int, skel_tree_t*>>* forest = il->computeSkeletons(firstLayer, SMfield, smThreshold);
    auto end = std::chrono::steady_clock::now();
    auto diff_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
    auto diff_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Computing skeletons took = " << diff_seconds << " s (" << diff_ms << " ms)" << std::endl;

    delete il;
    return forest;
}

string colorspace_to_string() {
    switch (c_space) {
    case GRAY:
        return string("Gray");
    case RGB:
        return string("RGB");
    case HSV:
        return string("HSV");
    case YCC:
        return string("YCbCr");
    default:
        return string("NONE! (FIX THIS)");
    }
}


float get_min_elem(FIELD<float>* im) {
    float min_elem = 1e5;
    int nPix = im->dimX() * im->dimY();
    float *c = im->data();
    float *end = im->data() + nPix;
    while (c != end)
        min_elem = min(min_elem, *c++);
    return min_elem;
}

void execute_gray_pipeline(FIELD<float>* im) {
    int clear_color = get_min_elem(im);
    int firstLayer = 0;
    vector<std::pair<int, skel_tree_t*>>* forest = execute_skeleton_pipeline(im, &firstLayer);
    PRINT(MSG_NORMAL, "Creating ImageWriter object...\n");
    ImageWriter iw(OUTPUT_FILE.c_str());
    if (c_space == COLORSPACE::NONE)
        c_space = COLORSPACE::GRAY;
    PRINT(MSG_NORMAL, "Using colorspace %s\n", colorspace_to_string().c_str());
    iw.writeHeader(im->dimX(), im->dimY(), c_space, clear_color,firstLayer,0,0);
    iw.write_image(forest);
    delete im;
}

void execute_color_pipeline(IMAGE<float>* im) {
    for (float* r = im->r.data(), *g = im->g.data(), *b = im->b.data(), *rend = im->r.data() + im->dimX() * im->dimY(); r < rend; ++r, ++g, ++b) {
        unsigned char r_prime = static_cast<unsigned char>((*r) * 255.0);
        unsigned char g_prime = static_cast<unsigned char>((*g) * 255.0);
        unsigned char b_prime = static_cast<unsigned char>((*b) * 255.0);
        unsigned char Y, Cb, Cr;
        float H, S, V;
        double max_val, min_val, delta;
        switch (c_space) {
        case COLORSPACE::YCC:
            Y  = min(max(0, round( 0.299  * r_prime + 0.587  * g_prime + 0.114  * b_prime      )), 255);
            Cb = min(max(0, round(-0.1687 * r_prime - 0.3313 * g_prime + 0.5    * b_prime + 128)), 255);
            Cr = min(max(0, round( 0.5    * r_prime - 0.4187 * g_prime - 0.0813 * b_prime + 128)), 255);
            r_prime = Y;
            g_prime = Cb;
            b_prime = Cr;
            break;
        case COLORSPACE::HSV:
            min_val = fmin(*r, fmin(*g, *b));
            max_val = fmax(*r, fmax(*g, *b));
            delta = max_val - min_val;
            H = 0.0;
            V = max_val;
            S = max_val > 1e-6 ? delta / max_val : 0.0f;

            if (S > 0.0f) {
                if (*r == max_val)       H = (0.0f + (*g - *b) / delta) * 60.0;
                else if (*g == max_val)  H = (2.0f + (*b - *r) / delta) * 60.0;
                else                H = (4.0f + (*r - *g) / delta) * 60.0;


                if (H < 0.0f)   H += 360.0f;
            }
            H = fmod(H, 360.0f) / 360.0f;
            r_prime = round(H * 255.0);
            g_prime = round(S * 255.0);
            b_prime = round(V * 255.0);
            break;
        default:
            // Do nothing and encode RGB
            break;
        }
        *r = r_prime, *g = g_prime, *b = b_prime;
    }

    FIELD<float> red_channel   = im->r;
    FIELD<float> green_channel = im->g;
    FIELD<float> blue_channel  = im->b;
    red_channel.writePGM("Y.pgm");
    green_channel.writePGM("CR.pgm");
    blue_channel.writePGM("CB.pgm");
    int num_layers_old;
    if (c_space == COLORSPACE::NONE)
        c_space = COLORSPACE::RGB;
    float r_min = get_min_elem(&(im->r));
    float g_min = get_min_elem(&(im->g));
    float b_min = get_min_elem(&(im->b));
    int Rfirst = 0;
    int Gfirst = 0;
    int Bfirst = 0;
    auto red_forest = execute_skeleton_pipeline(&red_channel, &Rfirst);
    time1 = 0;
    auto green_forest = execute_skeleton_pipeline(&green_channel, &Gfirst);
    time1 = 0;
    auto blue_forest = execute_skeleton_pipeline(&blue_channel, &Bfirst);
    PRINT(MSG_NORMAL, "Creating ImageWriter object...\n");
    ImageWriter iw(OUTPUT_FILE.c_str());
    PRINT(MSG_NORMAL, "Using colorspace %s\n", colorspace_to_string().c_str());
    iw.writeHeader(im->dimX(), im->dimY(), c_space, min(min(r_min, g_min), b_min), Rfirst, Gfirst, Bfirst);
    iw.write_color_image(red_forest, green_forest, blue_forest);
}

int main(int argc, char **argv) {
    // Increase the stacksize because default is not enough for some reason.
    const rlim_t kStackSize = 128 * 1024 * 1024;   // min stack size = 128 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0) {
        if (rl.rlim_cur < kStackSize) {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0) {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            } else {
                PRINT(MSG_NORMAL, "Set stack limit to %d\n", kStackSize);
            }
        }
    }

    /* Parse options */
    if (argc != 2) {
        fstream fhelp(". / doc / USAGE", fstream::in);
        cout << fhelp.rdbuf();
        fhelp.close();

        exit(-1);
    }
    Config config = Config(argv[1], NULL);
    setConfig(config);
    if (argc == 3) {
        OUTPUT_FILE = argv[2];
    }
    const char* filename = filename_stdstr.c_str();
    const char* saliencyMap = filename_sm.c_str();//wang

    PRINT(MSG_NORMAL, "Reading input file : % s\n", filename);
    
    // TODO(maarten): Rely on smt better than file type detection or expand the list
    bool is_color_image = (FIELD<float>::fileType((char*)filename) ==
                           FIELD<float>::FILE_TYPE::PPM);
    if (is_color_image) {
        IMAGE<float>* f = IMAGE<float>::read(filename); /* Read scalar field input */
        if (!f) {
            PRINT(MSG_ERROR, "Failed to read file.\n");
            exit(EXIT_FAILURE);
        }
        PRINT(MSG_NORMAL, "Executing color pipeline\n");
        //wang.
        SMfield = FIELD<float>::read(saliencyMap);

        if (!SMfield) {
            PRINT(MSG_ERROR, "Failed to read file.\n");
            exit(EXIT_FAILURE);
        }
        execute_color_pipeline(f);
    } else {
        PRINT(MSG_NORMAL, "Executing grayscale pipeline\n");
        FIELD<float>* field = FIELD<float>::read(filename);

        if (!field) {
            PRINT(MSG_ERROR, "Failed to read file.\n");
            exit(EXIT_FAILURE);
        }
         //wang
         SMfield = FIELD<float>::read(saliencyMap);

        if (!SMfield) {
            PRINT(MSG_ERROR, "Failed to read file.\n");
            exit(EXIT_FAILURE);
        }

        execute_gray_pipeline(field);
    }

    PRINT(MSG_NORMAL, "Done with everything.\n");
    return 0;
}




