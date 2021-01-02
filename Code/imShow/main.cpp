#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <unistd.h>
#include <tuple>
#include "include/main.hpp"
#include "include/io.hpp"
#include "include/messages.h"
#include "include/glslProgram.h"
#include "include/framebuffer.hpp"
#include "include/DataManager.hpp"
#include "include/ImageDecoder.hpp"
#include "lodepng/lodepng.h"
#include "include/image.h"
#include "fileio/fileio.hpp"

#define MIN(a,b) (a) < (b) ? (a):(b)
#define MAX(a,b) (a) > (b) ? (a):(b)
#define ARRAY_SIZE(arr) sizeof(arr)/sizeof(arr[0])


SHADER_TYPE SHADER = NORMAL;
//SHADER_TYPE SHADER = LOCAL_ALPHA;
// Possible values: GLOBAL_ALPHA, LOCAL_ALPHA, NORMAL;
/* Global alpha uses a distance transform on the entire layer.
 * Local alpha uses a distance for each disk (TESTING)
 * Normal does not perform additional interpolation.
 */

/* Output level: MSG_QUIET, MSG_ERROR, MSG_NORMAL, MSG_VERBOSE */
char * output_num; bool SaveMore = 0;//wang
int MSG_LEVEL = MSG_VERBOSE;
int WWIDTH = 0, WHEIGHT = 0;
int clear_color;
unsigned char* new_pixel_data = 0;
DataManager *dm, *dm_g, *dm_b;
bool is_color_image = false;
COLORSPACE c_space = COLORSPACE::NONE;
GLSLProgram *sh_render;
bool interpolate = true;
int display1;
using namespace std;
int firstLayerNum; int firstLayerRGB = 0; int time1;//wang
char *outFile = 0; /* Press 's' to make a screenshot to location specified by outFile. */

void saveOutput() {
    unsigned char *sdata = (unsigned char *) malloc(WWIDTH * WHEIGHT * 4);
    glReadPixels(0, 0, WWIDTH, WHEIGHT, GL_RGB, GL_UNSIGNED_BYTE, sdata);
    Texture2D *t = new Texture2D(WWIDTH, WHEIGHT, GL_RGB, GL_RGB, GL_UNSIGNED_BYTE);
    t->setData(sdata);
    t->saveAsPNG(outFile);
    free(sdata);
    delete t;
}

void initDataManager(const char *file) {
    ImageDecoder id;
    c_space = id.load(file);
    //PRINT(MSG_NORMAL, "11\n");
    is_color_image = (c_space != COLORSPACE::NONE && c_space != COLORSPACE::GRAY);
    
    PRINT(MSG_NORMAL, "Is color image? %s\n", is_color_image ? "Yes" : "No");
    dm = new DataManager();
    dm->setWidth(id.width);
    dm->setHeight(id.height);
    dm->initCUDA();
    dm->setData(id.getImage());
    clear_color = id.clear_color;
    dm->set_clear_color(id.clear_color);
    dm->set_gray_levels(id.get_gray_levels());
    dm->setFirstLayer(id.Rfirst);//wang

    if (is_color_image) {
        dm_g = new DataManager();
        dm_g->setWidth(id.width);
        dm_g->setHeight(id.height);
        dm_g->setData(id.getGChannel());
        dm_g->set_clear_color(id.clear_color);
        dm_g->set_gray_levels(id.get_g_levels());
        dm_g->setFirstLayer(id.Gfirst);//wang

        dm_b = new DataManager();
        dm_b->setWidth(id.width);
        dm_b->setHeight(id.height);
        dm_b->setData(id.getBChannel());
        dm_b->set_clear_color(id.clear_color);
        dm_b->set_gray_levels(id.get_b_levels());
        dm_b->setFirstLayer(id.Bfirst);//wang

    }
}

void initShader() {
    string runifs[] = {"alpha",
                       "layer"
                      };
    sh_render = new GLSLProgram("glsl/render.vert", "glsl/render.frag");
    sh_render->compileAttachLink();
    sh_render->bindUniforms(ARRAY_SIZE(runifs), runifs);
}

void draw_image(DataManager* a_dm) {
    int prev_layer = 0;
    int prev_layer4 = 0;//wang
    Texture2D* alpha = 0;
    glClear(GL_COLOR_BUFFER_BIT);
    vector<int> gray_levels = a_dm->get_gray_levels();
    int max_level = *std::max_element(gray_levels.begin(), gray_levels.begin() + gray_levels.size());


     //Draw the first white layer before other layers.
     firstLayerRGB ++;
     if (firstLayerRGB == 1) firstLayerNum = a_dm->getFirstLayer();
     else if (firstLayerRGB == 2) firstLayerNum = a_dm->getFirstLayer();
          else firstLayerNum = a_dm->getFirstLayer();

        alpha = a_dm->getWhiteMap();
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        sh_render->bind();
        glActiveTexture(GL_TEXTURE0);
        alpha->bind();

        glUniform1f(sh_render->uniform("layer"), firstLayerNum / 255.0);
        glBegin(GL_QUADS);
        glTexCoord2f(0.0, 1.0);
        glVertex2f(0, 0);
        glTexCoord2f(0.0, 0.0);
        glVertex2f(0, WHEIGHT);
        glTexCoord2f(1.0, 0.0);
        glVertex2f(WWIDTH, WHEIGHT);
        glTexCoord2f(1.0, 1.0);
        glVertex2f(WWIDTH, 0);
        glEnd();

        sh_render->unbind();
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);

    for (int inty : gray_levels) {
        time1++;
        //cout<<time1<<" time- level "<<inty<<endl;
        //if(inty == 255) inty = 254;  
        // PRINT(MSG_NORMAL, "%d -> %d/%d\n", prev_layer, inty, max_level);
        if (interpolate)
        {
            //Adjacent layer interpolation
            //alpha = a_dm->get_interp_layer(inty, prev_layer); prev_layer = inty;
            
            //Non-adjacent layer interpolation
            if(time1 % 2 == 0)
            {
                //cout<<time1<<" time- level "<<inty<<endl;
                alpha = a_dm->get_interp_layer(inty, prev_layer4);
                prev_layer4 = inty;
            }
            else
            {
                alpha = a_dm->get_interp_layer(inty, prev_layer);
                prev_layer = inty;
            }
            /**/

        }
        else
            alpha = a_dm->getAlphaMapOfLayer(inty);
            
            //string s = "reconstruct" + std::to_string(inty) + ".png";
                //alpha -> saveAsPNG(s);
           
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        /* SECOND PASS: Render using alpha map */
        sh_render->bind();
        glActiveTexture(GL_TEXTURE0);

        alpha->bind();
        glUniform1f(sh_render->uniform("layer"), inty / 255.0);
        glBegin(GL_QUADS);
        glTexCoord2f(0.0, 1.0);
        glVertex2f(0, 0);
        glTexCoord2f(0.0, 0.0);
        glVertex2f(0, WHEIGHT);
        glTexCoord2f(1.0, 0.0);
        glVertex2f(WWIDTH, WHEIGHT);
        glTexCoord2f(1.0, 1.0);
        glVertex2f(WWIDTH, 0);
        glEnd();
        sh_render->unbind();
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        //prev_layer = inty;
    }
}

// Convert some colorspace back to RGB
tuple<unsigned char, unsigned char, unsigned char> convert_color_space(unsigned char x, unsigned char y, unsigned char z) {
    int r, g, b;
    float r_prime, g_prime, b_prime;


    if (c_space == COLORSPACE::YCC) {
        r = min(max(0, round(x                      + 1.402  * (z - 128))), 255);
        g = min(max(0, round(x - 0.3441 * (y - 128) - 0.7141 * (z - 128))), 255);
        b = min(max(0, round(x + 1.772  * (y - 128)                     )), 255);
    } else if (c_space == COLORSPACE::HSV) {
        float h = (x / 255.0) * 360.0;
        float s = y / 255.0;
        float v = z / 255.0;
        // cout << h << " - " << s << " - " << v << endl;
        if (s == 0.0) {
            r_prime = g_prime = b_prime = z;
        } else {
            int hi = static_cast<int>(h / 60.0f);
            double f = h / 60.0f - hi;
            double p = v * (1.0f - s);
            double q = v * (1.0f - f * s);
            double t = v * (1.0f - (1.0f - f) * s);

            switch (hi) {
            case 0: r_prime = v; g_prime = t; b_prime = p; break;
            case 1: r_prime = q; g_prime = v; b_prime = p; break;
            case 2: r_prime = p; g_prime = v; b_prime = t; break;
            case 3: r_prime = p; g_prime = q; b_prime = v; break;
            case 4: r_prime = t; g_prime = p; b_prime = v; break;
            default: r_prime = v; g_prime = p; b_prime = q; break;
            }
        }
        r = round(r_prime * 255.0);
        g = round(g_prime * 255.0);
        b = round(b_prime * 255.0);
    } else {
        // RGB
        r = x, g = y, b = z;
    }
    return make_tuple(r, g, b);
}

void draw_color_image(DataManager* dm_r, DataManager* dm_g, DataManager* dm_b) {
    time1=1;//In order to be consistent with the time in imConvert file;
    draw_image(dm_r);

    FIELD<float>* chan1 = dm_r->get_texture_data();
    chan1->writePGM("r.pgm");

    time1=1;//In order to be consistent with the time in imConvert file;
    draw_image(dm_g);
    FIELD<float>* chan2 = dm_g->get_texture_data();
    chan2->writePGM("g.pgm");

    time1=1;//In order to be consistent with the time in imConvert file;
    draw_image(dm_b);
    FIELD<float>* chan3 = dm_b->get_texture_data();
    chan3->writePGM("b.pgm");


    new_pixel_data = (unsigned char*)(malloc(WWIDTH * WHEIGHT * 4 * sizeof(unsigned char)));
    for (int i = 0; i < WWIDTH; ++i) {
        for (int j = 0; j < WHEIGHT; ++j) {
            float a = chan1->value(i, j);
            float b = chan2->value(i, j);
            float c = chan3->value(i, j);
            auto rgb_tuple = convert_color_space(a, b, c);

            int id = j * WWIDTH + i;
            new_pixel_data[4 * id] = get<0>(rgb_tuple);
            new_pixel_data[4 * id + 1] = get<1>(rgb_tuple);
            new_pixel_data[4 * id + 2] = get<2>(rgb_tuple);
            new_pixel_data[4 * i + 3] = 255;
        }
    }
    glTexSubImage2D(GL_TEXTURE_2D, 0 , 0, 0, WWIDTH, WHEIGHT, GL_RGBA,
                    GL_UNSIGNED_BYTE, reinterpret_cast<GLvoid*>(new_pixel_data));

     glBegin(GL_QUADS);
    glTexCoord2d(0.0, 0.0);      glVertex2d(0.0,    WHEIGHT);
    glTexCoord2d(1.0, 0.0);      glVertex2d(WWIDTH, WHEIGHT);
    glTexCoord2d(1.0, 1.0);      glVertex2d(WWIDTH, 0.0);
    glTexCoord2d(0.0, 1.0);      glVertex2d(0.0,    0.0);
    glEnd();
}


void display(void) {
    dm->setClearColor();
    if (is_color_image) {
        dm_b->setClearColor();
        dm_g->setClearColor();
    }

    cout << is_color_image << endl;
    if (is_color_image) {
        draw_color_image(dm, dm_g, dm_b);
    } else {
       draw_image(dm);
    }
    //glFlush();  // Finish rendering
   // glutSwapBuffers();
    if (SaveMore) //wang.
    {
        stringstream ss;
        ss<<"output"<<output_num<<".png";
        outFile = const_cast<char*>(ss.str().c_str());
    }
    else  outFile = const_cast<char*>("output.png");

    saveOutput();
    glutDestroyWindow(display1);
}

void reshape(int w, int h) {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, WWIDTH, WHEIGHT, 0);
    glViewport(0, 0, WWIDTH, WHEIGHT);
    glMatrixMode(GL_MODELVIEW);
    return;
}

void idle(void) {
    usleep(1000);
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 's':
        outFile = const_cast<char*>("output.png");
        saveOutput();
        return;
        break;
    case 'q':
    case 27:
        exit(0);
        break;
    case 32:
        // TODO(maarten): interpolated image caching?
        interpolate = !interpolate;
        break;
    default:
        printf("Unsupported input: %c", key);
        fflush(stdout);
        break;
    }
    glutPostRedisplay();
}

int main(int argc, char *argv[]) {

    if (argc==3) //need to save lots of output images.
    {
        SaveMore = 1;
        output_num = argv[2];
    }

    /* Read meta data and set variables accordingly */
    initDataManager(argv[1]);
    WHEIGHT = dm->getHeight();
    WWIDTH = dm->getWidth();
    PRINT(MSG_VERBOSE, "Image dimensions: %d x %d\n", WWIDTH, WHEIGHT);


    /* Initialize GLUT */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(WWIDTH, WHEIGHT);
    display1 = glutCreateWindow(WINDOW_TITLE);

    /* Initialize GLEW */
    GLenum err = glewInit();
    if (GLEW_OK != err) {
        /* Problem: glewInit failed, something is seriously wrong. */
        PRINT(MSG_ERROR, "Error: %s\n", glewGetErrorString(err));
        exit(-1);
    }
    PRINT(MSG_NORMAL, "GLEW: Using GLEW %s\n", glewGetString(GLEW_VERSION));

    /* Set OPENGL states */
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    /* Set clear color one below the first layer that has points*/

    /* Set texture parameters */
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);


    glutDisplayFunc(display);
     
   // glutKeyboardFunc(keyboard);
    glutReshapeFunc(reshape);
   // glutIdleFunc(idle);


    initShader();
    dm->initFBO();
    dm->setAlphaShader(SHADER);
    if (is_color_image) {
        dm_g->initFBO();
        dm_g->setAlphaShader(SHADER);
        dm_b->initFBO();
        dm_b->setAlphaShader(SHADER);
    }

    
    glutMainLoop();
    return 0;
}