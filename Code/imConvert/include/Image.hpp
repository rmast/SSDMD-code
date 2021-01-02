#ifndef IMAGE_GUARD
#define IMAGE_GUARD

#include "CUDASkel2D/include/field.h"
#include "ImageWriter.hpp"
#include <set>
#include <vector>

class Image {
  private:
    /** VARIABLES **/
    float      islandThreshold;
    float             layerThreshold;
    int               DistinguishableInterval;
    double            *importance;
    //double            *UpperLevelSet;
    int               numLayers;
    int               nPix; /* Short for (dimX * dimY) */
    string            compress_method;
    int*  graylevels = nullptr;
    skel_tree_t *traceLayer(FIELD<float> *skel, FIELD<float> *dt);
    skel_tree_t *tracePath(int x, int y, FIELD<float> *skel, FIELD<float> *dt);
    coord2D_list_t *neighbours(int x, int y, FIELD<float> *skel);

  public:
    /** VARIABLES **/
    FIELD<float>   *im;

    /** CONSTRUCTORS **/
    Image(FIELD<float> *in);
    Image(FIELD<float> *in, float islandThresh, float importanceThresh, int GrayInterval);

    /** DESTRUCTOR **/
    ~Image();

    /** FUNCTIONS **/
    static void removeIslands(FIELD<float>*layer, unsigned int iThresh);
    void removeIslands(FIELD<float>* sm, float flagThreshold);
    void removeLayers();
    void calculateImportance();
    std::vector<std::pair<int, skel_tree_t*>>* computeSkeletons(int *firstLayer, FIELD<float>* sm, float smThreshold);
    void computeCUDASkeletons();
    void removeDuplicatePoints(FIELD<float> *imPrev, FIELD<float> *skP, FIELD<float> *imCur, FIELD<float> *skC);
    pair<int, int> find_closest_point(int i, int j, FIELD<float>* skelPrev);
    void bundle(FIELD<float>* skelCurr, FIELD<float>* skelPrev, FIELD<float>* currDT, FIELD<float>* prevDT, short* prev_skel_ft, int fboSize);
    double overlap_prune(FIELD<float>* skelPrev, FIELD<float>* currDT, FIELD<float>* prevDT);
};

#endif
