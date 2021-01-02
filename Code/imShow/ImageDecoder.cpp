/*
 * File:   ImageDecoder.cpp
 * Author: yuri
 *
 * Created on June 13, 2011, 1:32 PM
 */

#include <vector>
#include "./include/io.hpp"
#include "./include/ImageDecoder.hpp"
#include "./include/messages.h"
#include "./include/DataManager.hpp"
#include "../shared/FastAC/arithmetic_codec.h"
#include "squash/squash.h"
#include "fileio/fileio.hpp"
#include <assert.h>
#include <fstream>
#include <bitset>
#include <climits>


using namespace std;

ImageDecoder::ImageDecoder() {
    width = -1;
    height = -1;
}

ImageDecoder::ImageDecoder(const ImageDecoder& orig) {
}

int MAXR = 0;

ImageDecoder::~ImageDecoder() {
}

/* Add a point to the current datatype. This is used for both starting points, and neighbouring points.
 * Identify a starting point by setting "point" to 0. This cannot happen for neighbouring points, as that
 * implicates the radius is 0.*/
void ImageDecoder::addPoint(char point, int &x, int &y, int &r, path_t *path) {
    /* If we are a neighbouring point: decode the character and walk in the direction. Modify the original
     * x,y,r values. */
    // Layout of point is 0xxyyrrr where xx are the bits denoting delta_x [0, 1, 2]
    // yy denotes delta_y [0, 1, 2] and rrr denote delta_radius [0, 1, 2, 3, 4]
    int8_t dx = ((point >> 5) & 0x3) - 1;
    int8_t dy = ((point >> 3) & 0x3) - 1;
    int8_t dr = (point & 0x7) - 2;
    x += dx;
    y += dy;
    r += dr;


    // if (r < 0) { cerr << "Invalid radius size, must be >0 [r=" << r << "]" << endl; exit(-1);}

    /* Add to datatype */
    path->push_back(coord3D_t(x, y, r));

}

const char* ImageDecoder::get_comp_method_string(COMPRESS_MODE mode) {
    switch (mode) {
    case COMPRESS_ZLIB:
        return "zlib";
        break;
    case COMPRESS_LZHAM:
        return "lzham";
        break;
    case COMPRESS_BSC:
        return "bsc";
        break;
    case COMPRESS_CSC:
        return "csc";
        break;
    case COMPRESS_LZMA:
        return "xz";
        break;
    case COMPRESS_LZMA2:
        return "lzma2";
        break;
    case COMPRESS_BROTLI:
        return "brotli";
        break;
    case COMPRESS_ZPAQ:
        return "zpaq";
        break;
    case COMPRESS_BZIP2:
        return "bzip2";
        break;
    case COMPRESS_ZSTD:
        return "zstd";
        break;
    default:
        // Uses some unsupported method or file is mangled.
        return NULL;
        break;
    }
}

/* Decode LZMA. */
int ImageDecoder::decompress_data(unsigned char *compr_data, int length, unsigned char **out, unsigned int *outlength) {
    unsigned char *dat_it = compr_data;
    unsigned char *decompr_data;
    size_t decompr_length = -1;
    size_t compr_length = length;
    /* First byte stores the compression method. */
    int compression_method = (int)dat_it[0];
    const char* comp_method_name = get_comp_method_string((ImageDecoder::COMPRESS_MODE)compression_method);
    dat_it += 1; compr_length -= 1;

    /* First 4 bytes store the size of the uncompressed data. */
    decompr_length = *((unsigned int *) dat_it);
    decompr_data = new unsigned char[decompr_length];

    dat_it += 4; compr_length -= 4;

    // TODO: FIX THIS SO IT ALWAYS WORKS!
    squash_set_default_search_path("../shared/squash/plugins");

    SquashCodec* codec = squash_get_codec(comp_method_name);
    if (codec == NULL) {
        fprintf (stderr, "Unable to find algorithm '%s'.\n", comp_method_name);
        delete [] decompr_data;
        return EXIT_FAILURE;
    }


    PRINT(MSG_NORMAL, "%s Input size: %f KB, Uncompressed output: %fKB\n", comp_method_name, length / 1024.0, decompr_length / 1024.0);
    SquashStatus res = squash_codec_decompress (codec,
                       &decompr_length, (uint8_t*) decompr_data,
                       compr_length, dat_it, NULL);

    switch (res) {
    case SQUASH_OK:
        cout << comp_method_name << " decompression went OK." << endl;
        break;
    // TODO(maarten): MORE CASES HERE!
    default:
        cout << comp_method_name << " returned unknown return value" << endl;
    }
    /* Fill return values */
    *out = decompr_data;
    *outlength = decompr_length;

    free (compr_data);
    return res;
}

bool ImageDecoder::decode_layer(image_t** image_ref, int intensity) {
    // uint8_t x1 = READUINT8(dat_it);
    // uint8_t x2 = READUINT8(dat_it);
   if (intensity == 255) intensity = 254;//Wang. I don't know why when the intensity equals to 255, it doesn't work, so here I just change it to 254 since it isn't distinguishable.
    image_t* image = *image_ref;
    path_t* layer = (*image)[intensity];
    // uint16_t numPaths =  (x1 << 8) + x2;
    uint16_t numPaths = READUINT16(dat_it);
    uint8_t current;
    if (numPaths == 0xFFFF)
        return false;

    of_uncompressed << "Intensity " << +intensity << " - " << "Num children " << +numPaths << endl;
    
    for (int i = 0; i < numPaths; ++i) {
     
        path_t* path = new path_t();
        /* Read first -full- point */
        /* Seriously, y tho*/
        uint8_t x1 = READUINT8(dat_it);
        uint8_t x2 = READUINT8(dat_it);
        int x = (x1 << 8) + x2;
        uint8_t y1 = READUINT8(dat_it);
        uint8_t y2 = READUINT8(dat_it);
        int y = (y1 << 8) + y2;
        uint8_t r1 = READUINT8(dat_it);
        uint8_t r2 = READUINT8(dat_it);
        int r = (r1 << 8) + r2;
        path->push_back(coord3D_t(x, y, r));
      
        bool end = false;
        while (!end) {
            of_uncompressed << x << " - " << y << " - " << r << endl;
            current = READUINT8(dat_it);
            if (current == END_TAG) {
                of_uncompressed << "End" << endl;
                end = true;
            }

            /* End of branch? */
            else if (current == FORK_TAG) {
                uint8_t goBack_1 = READUINT8(dat_it);
                uint8_t goBack_2 = READUINT8(dat_it);
                uint16_t goBack = (goBack_1 << 8) + goBack_2; // TODO(maarten): why is this necessary?
                of_uncompressed << "Fork - " << (goBack) << endl;
                for (unsigned int q = 0; q < goBack; ++q) {
                    layer->push_back((path->back()));
                    path->pop_back();
                }
                x = path->back().first;
                y = path->back().second;
                r = path->back().third;
            } else {
                if ((current < 128)) {
                    addPoint(current, x, y, r, path);
                } else {
                    if (current == 128) {
                        if (bits_dx == 0) {
                            int8_t dx_8 = READINT8(dat_it);
                            x += dx_8;
                        } else {
                            uint8_t dx_1_16 = READUINT8(dat_it);
                            uint8_t dx_2_16 = READUINT8(dat_it);
                            x += ((dx_1_16 << 8) + dx_2_16);
                            if (x > width)
                                x -= (1 << 16);
                        }
                        if (bits_dy == 0) {
                            int8_t dy_8 = READINT8(dat_it);
                            y += dy_8;
                        } else {
                            uint8_t dy_1_16 = READUINT8(dat_it);
                            uint8_t dy_2_16 = READUINT8(dat_it);
                            y += ((dy_1_16 << 8) + dy_2_16);
                            if (y > width)
                                y -= (1 << 16);

                        }
                        if (bits_dr == 0) {
                            int8_t dr_8 = READINT8(dat_it);
                            r += dr_8;
                        } else {
                            uint8_t dr_1_16 = READUINT8(dat_it);
                            uint8_t dr_2_16 = READUINT8(dat_it);
                            r += ((dr_1_16 << 8) + dr_2_16);
                            if (r > width)
                                r -= (1 << 16);

                        }
                    } else {
                        uint16_t num = (current << 8) + READUINT8(dat_it);
                        int8_t dr = ((num & 0x1F)) - 15;
                        int8_t dy = (((num >> 5) & 0x1F)) - 15;
                        int8_t dx = (((num >> 10) & 0x1F)) - 15;
                        x += dx;
                        y += dy;
                        r += dr;
                    }

                    if (r < 0) { cerr << "Invalid radius size, must be >0 [r=" << r << "]" << endl; exit(-1);}

                    /* Add to datatype */
                    path->push_back(coord3D_t(x, y, r));
                }
            }

        }
        
        /* Copy path to the layer of the image */
        for (unsigned int i = 0; i < path->size(); ++i) {
           /* if (intensity ==249)
            {
                if(i == 0){     //just for debug
                PRINT(MSG_NORMAL, "pathx %d\n", (*path)[1].first);
                PRINT(MSG_NORMAL, "pathy %d\n", (*path)[1].second);
                PRINT(MSG_NORMAL, "pathz %d\n", (*path)[1].third);
                }
            } */
                
            layer->push_back((*path)[i]);
        }

    }
    //PRINT(MSG_NORMAL, "33\n");
    return true;
}

image_t* ImageDecoder::decode_plane(vector<int>& levels) {
    image_t* img = new image_t();
   
    levels.clear();
    int num_levels = READUINT8(dat_it);
    for (int i = 0; i < num_levels; ++i) {
        int level = READUINT8(dat_it);
        if (i == 0)
            levels.push_back(level);
        else
            levels.push_back(levels.at(i - 1) + level);
    }
    for (int i = 0; i < 0xff; ++i) {
        path_t *p = new path_t();
        img->push_back(p);
    }
    of_uncompressed.open("uncompressed.sir", ios_base::out | ios::binary);
    //PRINT(MSG_NORMAL, "1\n");
    for (auto intensity : levels) {
        decode_layer(&img, intensity);
    }
    of_uncompressed.close();
    return img;
}


/* Load SIR file, decode compressed data, and read all disks. */
COLORSPACE ImageDecoder::load(const char *fname) {
    unsigned int nEl = readFile<unsigned char>(fname, &data);
    if (nEl == 0) {
        PRINT(MSG_ERROR, "Could not open file for Image Decoder.\n");
        exit(1);
    }

    /* Decode compressed data, and overwrite old data with the decompressed data*/
    decompress_data(data, nEl, &data, &nEl);

    /* Initialize iterator */
    dat_it = data;

    int version = READUINT16(dat_it);
    if (version != READER_FILE_VERSION_NUMBER) {
        cerr << "Incorrect version. Can read: " <<  READER_FILE_VERSION_NUMBER << " , encountered: " <<  version << endl;
    }

    for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
            for (int r = -2; r <= 2; ++r)
                point_list.push_back(make_tuple(x, y, r));


    /* Get Width and Height */
    width =  READUINT16(dat_it);
    height = READUINT16(dat_it);
    COLORSPACE colorspace = (COLORSPACE)READUINT8(dat_it);
    clear_color = READUINT8(dat_it);
    cout << "Colorspace: " << colorspace << endl;
    if (colorspace == COLORSPACE::GRAY) {
        Rfirst = READUINT8(dat_it);
        cout << "no colors" << endl;
        r_image = decode_plane(gray_levels);
        
    } else {
        Rfirst = READUINT8(dat_it);
        Gfirst = READUINT8(dat_it);
        Bfirst = READUINT8(dat_it);
        cout << "colors" << endl;
        r_image = decode_plane(gray_levels);
        g_image = decode_plane(g_gray_levels);
        b_image = decode_plane(b_gray_levels);
    }
    return colorspace;
}
