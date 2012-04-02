/*
 ============================================================================
 Name        : main.c
 Author      : Tom Lin (lead programmer)
 Version     : CHECK PARAMS.H
 Copyright   : SORA-SYSTEMS
 Description : CHECK PARAMS.H
 ============================================================================
 */
#include "params.h"
#include "proto.h"
#include "image.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

//#include "lib_matrix.c"
#define PI 3.141592
#define OPT 1




int main(int argc, char *argv[]) {

	//YSH NEED TO BE INT

	int image_main_ysh[2][3][Y_SIZE][X_SIZE];
	unsigned char image_main[2][3][Y_SIZE][X_SIZE];
	unsigned char image_main_ybuf[2][3][Y_SIZE][X_SIZE];
	unsigned char image_buffer[2][Y_SIZE][X_SIZE];
	unsigned char image_marker[Y_SIZE][X_SIZE];
	unsigned char image_thresh[2][Y_SIZE][X_SIZE];
	//unsigned char image_label[2][Y_SIZE][X_SIZE];
	unsigned char image_feature[2][Y_SIZE][X_SIZE];
	unsigned char image_sub[2][3][Y_SIZE][X_SIZE];
	unsigned char image_out[3][Y_SIZE][X_SIZE];
	unsigned char image_test[2][3][Y_SIZE][X_SIZE];
	int image_test_ysh[2][3][Y_SIZE][X_SIZE];
    int c;
    int i;
/*
	int array_size = 16*16;
    t_array_index ***position_neighbor;
    *position_neighbor =  malloc(sizeof ***position_neighbor * array_size * 3 * array_size);
    //position_neighbor = malloc((array_size * sizeof(t_array_index)) + (3 * array_size * array_size * sizeof(t_array_index)) + (3 * array_size * sizeof(t_array_index)));


    for (c = 0; c < 3; c++)
    {
    	position_neighbor[c] =  malloc(sizeof(t_array_index) * array_size * array_size);
        for (i = 0; i < array_size; i++)
        {
        	position_neighbor[c][i] =  malloc(sizeof(t_array_index) * array_size);
        }
    }



    int neighbor_range = 8;
    double neighbor_coefficient_matrix[c][array_size][array_size][neighbor_range+1];

    t_array_index ***neighbor_coefficient_deviation;
    *neighbor_coefficient_deviation = malloc(sizeof ***neighbor_coefficient_deviation * array_size * 3 * array_size);
    for (c = 0; c < 3; c++)
    {
    	neighbor_coefficient_deviation[c] = malloc(sizeof(t_array_index) * array_size * array_size);
        for (i = 0; i < array_size; i++)
        {
        	neighbor_coefficient_deviation[c][i] = malloc(sizeof(t_array_index) * array_size);
        }
    }
*/

	//unsigned char image_general[2][3][Y_SIZE][X_SIZE];
	//char *image_path = "miyako.bmp";
	//char *image_0_path = "miyako_FFT.bmp";
	//char *image_1_path = "miyako_FFT_flt.bmp";
	//image_path = "MIYAKO.bmp";
	clock_t start = clock();
	//clock_t start_read = clock();
	//clock_t start_fft = clock();
/*
 * VAR DEFINE
 */
	//LISTING DEFAULT
	int threshold_rgb = 128; //90?
	int threshold_ysh = 128; //90?
	double threshold_ff = 0.3; // *100 // OPTIMAL 0.3?
	int lim = (Y_SIZE+X_SIZE)/2;
	int a = 0;  //OPTIMAL 0!
	int b = 32; //OPTIMAL 128?
	int clrv = 1; // 0 FOR ONE, VICE VERSA => CLRV = N+1
	//DEBUG ONLY

	//debug(str);


	char *path_fft_a = join(join(join(CORE_CODE_NAME, "_FFT_A_"), tochar(TIME)), ".bmp");
	char *path_fft_b = join(join(join(CORE_CODE_NAME, "_FFT_B_"), tochar(TIME)), ".bmp");
	char *path_fft_x = join(join(join(CORE_CODE_NAME, "_FFT_X_"), tochar(TIME)), ".bmp");
	char *path_fft_m = join(join(join(CORE_CODE_NAME, "_FFT_M_"), tochar(TIME)), ".bmp");

	char *path_ff_a = join(join(join(CORE_CODE_NAME, "_FF_A_"), tochar(TIME)), ".bmp");
	char *path_ff_b = join(join(join(CORE_CODE_NAME, "_FF_B_"), tochar(TIME)), ".bmp");
	char *path_ff_x = join(join(join(CORE_CODE_NAME, "_FF_X_"), tochar(TIME)), ".bmp");


	char *path_thr_a = join(join(join(CORE_CODE_NAME, "_THR_A_"), tochar(TIME)), ".bmp");
	char *path_thr_b = join(join(join(CORE_CODE_NAME, "_THR_B_"), tochar(TIME)), ".bmp");

	char *path_src_a = join(join(join(CORE_CODE_NAME, "_SRC_A_"), tochar(TIME)), ".bmp");
	char *path_src_b = join(join(join(CORE_CODE_NAME, "_SRC_B_"), tochar(TIME)), ".bmp");

	int debug = 1;

	//DEBUG ONLY END

/*
 * VAR DEFINE END
 */


/*
 * PARAM READ
 */


int ch;
opterr = 0;
/*
 * PARAM DEFINITION
 * s -> array size
 * m -> processing method
 * t -> RGB compare threshold
 * y -> YSH compare threshold
 * f -> feature threshold
 * c -> convolution amplifier
 * r -> FFT filter radius
 * d -> compared image, debug output only
 * e -> eclipse debug(read from file)
 * a -> debug image path a
 * b -> debug image path b
 * PARAM DEFINITION END
 */
while((ch=getopt(argc, argv, "s:m:t:y:c:r:d:e:va:b:")) != -1)
	{
	switch(ch)
	{
	case 's':
		X_SIZE = atol(optarg);
		Y_SIZE = atol(optarg);
		lim = atol(optarg);
		break;
	case 'v':
			//#define CORE_NAME	SORA_CORE
			//#define CORE_DEV_STATUS	PRE
			//#define CORE_VERSION	0.1.4a
			//#define CORE_CODE_NAME	KIRINO
			//#define CORE_PLATFORM	LINUX_UNIX
			printf("%s;%s;%s;%s;%s\n", CORE_NAME, CORE_DEV_STATUS,
					CORE_VERSION, CORE_CODE_NAME, CORE_PLATFORM);
			return EXIT_SUCCESS;
			break;
	case 'm':

		break;
	case 't':
		threshold_rgb = atol(optarg);
		break;
	case 'y':
		threshold_ysh = atol(optarg);
		break;
	case 'f':
		threshold_ff = atof(optarg);
		break;
	case 'c':
		clrv = atol(optarg);
		break;
	case 'r':
		b = atol(optarg);
		break;
	case 'd':
		path_fft_x = optarg;
		break;
	case 'e':

		debug = atol(optarg);
		break;
	case 'a':
		debug_a = optarg;
		//char *debug_b = "ko.bmp";
		//debug = atol(optarg);
		break;
	case 'b':
		debug_b = optarg;
		//char *debug_b = "ko.bmp";
		//debug = atol(optarg);
		break;
	default:

		break;
	}
	}
if(b >= (lim/2))
	b = lim/2;

/*
 * PARAM READ END
 */



/*
	int i;
	int j;
	//unsigned char *image_buf;
	for (i = 0; i < Y_SIZE; i++)
			for (j = 0; j < X_SIZE; j++) {
				// *(image_buf + 3*(X_SIZE*i + j)    ) = image_sub[2][Y_SIZE-i-1][j];
				// *(image_buf + 3*(X_SIZE*i + j) + 1) = image_sub[1][Y_SIZE-i-1][j];
				// *(image_buf + 3*(X_SIZE*i + j) + 2) = image_sub[0][Y_SIZE-i-1][j];
				//printf("%d:%d:%d|",image_sub[0][i][j],image_sub[1][i][j],image_sub[2][i][j]);
			}
*/

	/*
	int x;
	int y;
	for(x = 0; x <= X_SIZE; x++)
	{
		for(y = 0; y <= Y_SIZE; y++)
		{
			printf("::%d::",image_sub[0][y][x]);
			//main code here!
		}

	}
	*/


/*
 * DATA READ
 */
	//DECL
if(debug == 0)
{
/*
	char *delim3 = "|";
	char *delim2 = "#";
	char *delim1 = ":";
	char *rtn;
	char *ptn;
	char *btn;
	char *rest;
	char *xest;
	char *yest;
	int xi = 0;
	int yi = 0;
	int ii = 0;
	int ci = 0;
	long fslen = (	((lim*lim)*3)*3+ ((lim*lim)*3)*3	) *2+lim;
	//char c[fslen];

	char *c = (char *)malloc(fslen); // allocate memory
	//(unsigned char *)malloc((size_t)fslen);

	//fread(image_buf, sizeof(unsigned char), (size_t)(long)X_SIZE*Y_SIZE*3, fp);


	fgets(c, fslen, stdin);






	//DECL END
if(c != NULL)
{
	btn = strtok_r(c, delim3,&yest);
	while(btn != NULL)
	{
		ptn = strtok_r(btn, delim2,&xest);
		while(ptn != NULL)
			{
				rtn = strtok_r(ptn, delim1,&rest);
				while(rtn != NULL)
				{
					if(yi <= lim)
					{
						image_test[ii][ci][xi][yi] = atol(rtn);
						//printf("|X=%d, Y=%d ; split char x = %s|\n", xi, yi,rtn);
						xi++;
					}
					if(xi >= lim)
					{
						yi++;
						xi = 0;
					}
					rtn = strtok_r(NULL, delim1, &rest);
				}
				ptn = strtok_r(NULL, delim2, &xest);
				ci++;
				xi = 0;
				yi = 0;
			}
			btn = strtok_r(NULL, delim3, &yest);
			ii++;
			xi = 0;
			yi = 0;
			ci = 0;
		}
		xi = 0;
		yi = 0;
		ci = 0;
		ii = 0;
}
//free(btn);
//free(ptn);
//free(rtn);
free(c);
//free(yest);
//free(xest);
//free(rest);
//free(&c);

*/
	read_image(image_test, lim);
}
else if(debug  == 1)
{
	//image_test[ii][ci][xi][yi]
	read_bmp_color(image_test[0], debug_a); // <-
	read_bmp_color(image_test[1], debug_b); // <-
}


// 8 or 24 only!
ImageData* img = createImage(256,256,24);
//void disposeImage(ImageData *img);
readBMPfile("miyako.bmp",&img);
writeBMPfile("test_img.bmp",img);


/*
 * DATA READ END
 */
//-----------------------------------------//
assign(image_main[0][0], image_test[0][0]);// <-
assign(image_main[0][1], image_test[0][1]);// <-
assign(image_main[0][2], image_test[0][2]);// <-
//-----------------------------------------//
assign(image_main[1][0], image_test[1][0]);// <-
assign(image_main[1][1], image_test[1][1]);// <-
assign(image_main[1][2], image_test[1][2]);// <-
//-----------------------------------------//
rgb_to_ysh(image_main[0],image_main_ysh[0]);// ->
rgb_to_ysh(image_main[1],image_main_ysh[1]);// ->
//::
y_image(image_main_ysh[0][0], image_main_ybuf[0][0]);// ->
y_image(image_main_ysh[1][0], image_main_ybuf[1][0]);// ->
//::
//threshold_dynamic(image_main_ybuf[0][0],image_thresh[0], 1);// ->
//threshold_dynamic(image_main_ybuf[1][0],image_thresh[1], 1);// ->
threshold_discrim(image_main_ybuf[0][0],image_thresh[0], 1);// ->
threshold_discrim(image_main_ybuf[1][0],image_thresh[1], 1);// ->
//::
median(image_thresh[0],image_buffer[0]);// ->
median(image_thresh[1],image_buffer[1]);// ->
//::
assign(image_thresh[0], image_buffer[0]);// <-
assign(image_thresh[1], image_buffer[1]);// <-
//::
int cnt_a;
int cnt_b;
int cc;


//char *buf;
/*
char text_buf_a[65536];
char text_buf_b[65536];
char buf_f[65536];
unsigned char image_label_out_a[Y_SIZE][X_SIZE];
unsigned char image_label_out_b[Y_SIZE][X_SIZE];
double size_a[HIGH], ratio_a[HIGH], length_a[HIGH];
double size_b[HIGH], ratio_b[HIGH], length_b[HIGH];
double size_f[HIGH], ratio_f[HIGH], length_f[HIGH];
*/
/*
unsigned char image_ffl[2][Y_SIZE][X_SIZE];
cnt_a = labeling(image_thresh[0], image_label[0]); //<<
cnt_b = labeling(image_thresh[1], image_label[1]); //<<
assign(image_ffl[0], image_label[0]);// <-
assign(image_ffl[1], image_label[1]);// <-
*/
//::

//labeling(image_thresh[0], image_label[0], *cnt, *buf);
//features(image_label[0], image_feature[0], cnt_a, size_a, length_a, ratio_a, text_buf_a); // ->.
//features(image_label[1], image_feature[1], cnt_b, size_b, length_b, ratio_b, text_buf_b); // ->.

//printf("FEATURE A:%d\n", ;

//printf("CNT A:%d\n", cnt_a);
//printf("CNT B:%d\n", cnt_b);
//printf("\nA=>\n%s\n<=A\n", text_buf_a);
//printf("\nB=>\n%s\n<=B\n", text_buf_b);

/*
assign(image_general[0][0], image_test[0][0]);// <-
assign(image_general[0][1], image_test[0][1]);// <-
assign(image_general[0][2], image_test[0][2]);// <-
//-----------------------------------------//
assign(image_general[1][0], image_test[1][0]);// <-
assign(image_general[1][1], image_test[1][1]);// <-
assign(image_general[1][2], image_test[1][2]);// <-
*/

//print_image(image_general[0], 1);
double diff_seg[4][3];

rgb_to_ysh(image_main[0],image_test_ysh[0]);// ->
rgb_to_ysh(image_main[1],image_test_ysh[1]);// ->

int ysh = 0;
//segment; mode; border
if(ysh)
	general_segmentation(image_test_ysh[0], image_test_ysh[1], 8, 0, 7, diff_seg);
else
	general_segmentation(image_test[0], image_test[1], 8, 0, 0, diff_seg);

//print_image(image_general[0],0);
//double diff_general = general_distance(image_general[0], image_general[1]);

double rslt;
int label[HIGH];
/*
features_compare(
	image_ffl[0],
	image_ffl[1],
	image_label_out_a,
	image_label_out_b,
	cnt_a, cnt_b, size_f,
	length_f, ratio_f, buf_f, &rslt, label, &cc, threshold_ff); // ->.

label_masking(image_label[0], image_marker, label, cc); // ->
*/

//printf("\nCC=>\n%s\n<=CC\n", buf_f);
//double *crt_rtn = &rslt;
//printf("=>LAB: %d \n", &label);

//printf("%f;", rslt);


if(debug  == 1)
{
if(write_bmp_mono(image_marker, path_ff_x) != -1)
	{debug("FF_X OK!");}
else
{debug("FF_X FAIL!");}
if(write_bmp_mono(image_feature[0], path_ff_a) != -1)
	{debug("FF A OK!");}
else
{debug("FF FAIL!");}
if(write_bmp_mono(image_feature[1], path_ff_b) != -1)
	{debug("FF B OK!");}
else
{debug("FF FAIL!");}
if(write_bmp_mono(image_buffer[0], path_thr_a) != -1)
	{debug("THR A OK!");}
else
{debug("THR FAIL!");}
if(write_bmp_mono(image_buffer[1], path_thr_b) != -1)
	{debug("THR B OK!");}
else
{debug("THR FAIL!");}

write_bmp_color(image_test[0], path_src_a);
write_bmp_color(image_test[1], path_src_b);
}



//FINISH THE FEATURE DETECTION HERE.......




/*
 * FFT FILTER/FFT
 */
if (((fftfilter_ff(image_test[0][0], image_sub[0][0],a,b) != -1) &&
	(fftfilter_ff(image_test[0][1], image_sub[0][1],a,b) != -1) &&
	(fftfilter_ff(image_test[0][2], image_sub[0][2],a,b) != -1)) &&
	((fftfilter_ff(image_test[1][0], image_sub[1][0],a,b) != -1) &&
	(fftfilter_ff(image_test[1][1], image_sub[1][1],a,b) != -1) &&
	(fftfilter_ff(image_test[1][2], image_sub[1][2],a,b) != -1))
	)
	{
		//printf("FFT OK!\n");
	}
/*
 * FFT FILTER/FFT END
 */

/*
 * COMPARE
 */
unsigned char image_dist[2][3][Y_SIZE][X_SIZE];
unsigned char image_rslt[3][Y_SIZE][X_SIZE];
int image_ysh[2][3][Y_SIZE][X_SIZE];
int cl,j,fl;
int true_val=0,false_val=0,total_val=0;

if(debug  == 1)
{
write_bmp_color(image_sub[0], path_fft_a);
write_bmp_color(image_sub[1], path_fft_b);
}


unsigned char image_fftm[3][Y_SIZE][X_SIZE];
unsigned char image_fftr[3][Y_SIZE][X_SIZE];

image_multiplication_color(image_sub[0], image_sub[1], image_fftm);

/*
fftimage_ff(image_fftm[0],image_fftr[0], 0, 0);
fftimage_ff(image_fftm[1],image_fftr[1], 0, 0);
fftimage_ff(image_fftm[2],image_fftr[2], 0, 0);
*/

fftimage_ff(image_sub[1][0],image_fftr[0], 0, 0);
fftimage_ff(image_sub[1][1],image_fftr[1], 0, 0);
fftimage_ff(image_sub[1][2],image_fftr[2], 0, 0);

fftfeature_ff(image_test[0][0], image_test[1][0], image_fftr[0], a, b);
fftfeature_ff(image_test[0][1], image_test[1][1], image_fftr[1], a, b);
fftfeature_ff(image_test[0][2], image_test[1][2], image_fftr[2], a, b);


write_bmp_color(image_fftr, path_fft_m);


//if((WRITE_FFT0 != -1) && (WRITE_FFT1 != -1))
//{//check FFT
	//printf("FFT WRITE OK!\n");

	//printf("TOTAL Time elapsed: %fs\n", ((double)clock() - start) / CLOCKS_PER_SEC);}
	rgb_to_ysh(image_sub[0],image_ysh[0]);//ysh convert image 0 ->
	rgb_to_ysh(image_sub[1],image_ysh[1]);//ysh convert image 1 ->
	for (i = 0; i < Y_SIZE; i++)
	{//y coordinates
		for (j = 0; j < X_SIZE; j++)
		{//x coordinates
			/*
			 * EXECUTE COMPARE {
			 */

						//if ((long)sqrt((double)((j-X_SIZE/2)*(j-X_SIZE/2) + (i-Y_SIZE/2)*(i-Y_SIZE/2))) <= b)
			//if((image_sub[0][1][i][j]  != 0))
			if((long)sqrt((double)((j-X_SIZE/2)*(j-X_SIZE/2) + (i-Y_SIZE/2)*(i-Y_SIZE/2))) <= b)
			{// 0 CHECK START
				total_val++;
				for(cl = 0;cl <= 2; cl++)
				{
					//RGB layer
					//exclude 0(!filtered)
					//RGB DISTANCE
					image_dist[0][cl][i][j] = fabs(image_sub[0][cl][i][j]-image_sub[1][cl][i][j]);
					//YSH DISTANCE
					image_dist[1][cl][i][j] = fabs(image_ysh[0][cl][i][j]-image_ysh[1][cl][i][j]);
					//DISTANCE WRITE
				}
			//if((image_sub[0][1][i][j]  != 0)){
				if(
						/*(image_sub[0][1][i][j]  != 0) && */
						(

						((image_dist[0][0][i][j]
						+
						image_dist[0][1][i][j]
						+
						image_dist[0][2][i][j])
						<=
						threshold_rgb)

						||

						((image_dist[1][0][i][j] <= threshold_ysh)
						&&
						(image_dist[1][1][i][j] <= threshold_ysh)
						&&
						(image_dist[1][2][i][j] <= threshold_ysh))

						)
					)
				{
					//int fl;
					true_val++;
					for(fl = 0; fl <= 2; fl++)
					{
						image_rslt[fl][i][j] = 0;

					}
				}
				else
				{
					//int fl;
					false_val++;
					for(fl = 0; fl <= 2; fl++)
					{
						image_rslt[fl][i][j] = 255;
					}
				}
			}
			else
			{
				for(fl = 0; fl <= 2; fl++)
				{
					image_rslt[fl][i][j] = 255;
				}
			}// O CHECK END
					/*
					 * EXECUTE COMPARE END }
					 */
		}//x coordinate end
	}//y coordinate end
//}// check FFT end
/*
 * COMPARE END
 */
/*
double per_val = (double)true_val/total_val;
printf("******\n");
printf("PERCENTAGE SIMILAR: %f%%,\n", per_val*100);
printf("PIXEL TOTAL: %d,\n", total_val);
printf("PIXEL SIMILAR: %d,\n", true_val);
printf("PIXEL DIFFER: %d,\n", false_val);
printf("******\n");
*/
int clr;
//convolute(image_rslt[0], image_out[0], b);
//convolute(image_rslt[1], image_out[1], b);
//convolute(image_rslt[2], image_out[2], b);
//param = clrv
//loop value  = clr


/*
 * AMPLIFY CONVOLUTION
 */
clrv--;
unsigned char image_conv[clrv+1][3][Y_SIZE][X_SIZE];
assign(image_conv[0][0], image_rslt[0]); //<-
assign(image_conv[0][1], image_rslt[1]); //<-
assign(image_conv[0][2], image_rslt[2]); //<-
for(clr = 0;clr <= clrv;clr++)
{
	convolute(image_conv[clr][0], image_conv[clr+1][0], b);
	convolute(image_conv[clr][1], image_conv[clr+1][1], b);
	convolute(image_conv[clr][2], image_conv[clr+1][2], b);
}
assign(image_out[0], image_conv[clrv+1][0]); //<-
assign(image_out[1], image_conv[clrv+1][1]); //<-
assign(image_out[2], image_conv[clrv+1][2]); //<-
/*
 * AMPLIFY CONVOLUTION END
 */

/*
 * CALCULATE PX
 */

int c_true = 0;
int c_false = 0;
int c_total = 0;
for (i = 0; i < Y_SIZE; i++)
{//y coordinates
	for (j = 0; j < X_SIZE; j++)
	{//x coordinates
		if((long)sqrt((double)((j-X_SIZE/2)*(j-X_SIZE/2) + (i-Y_SIZE/2)*(i-Y_SIZE/2))) <= b)
		{
			c_total++;
			if((image_out[0][i][j] && image_out[1][i][j] && image_out[2][i][j]) == 0)
			{
				c_true++;
			}
			else
			{
				c_false++;
			}
		}
	}
}
/*
 * CALCULATE PX END
 */
//printf("|||||||\n");
//double diff_fft = general_segmentation(image_sub[0], image_sub[1],64);
double c_per = (double)c_true/c_total;
//printf("FFT: %f", c_per);
printf("%f;", c_per);
/*
printf("PERCENTAGE SIMILAR: %f%%,\n", c_per*100);
printf("PIXEL TOTAL: %d,\n", c_total);
printf("PIXEL SIMILAR: %d,\n", c_true);
printf("PIXEL DIFFER: %d,\n", c_false);
*/

//int c;
double deviation_unity;
for(c = 0;c < 3;c++)
{
	deviation_unity += diff_seg[1][c];
	//printf("D%d = %f\n",c,diff_seg[1][c]);
}

deviation_unity/=3;

double average_unity;
for(c = 0;c < 3;c++)
{
	average_unity += diff_seg[2][c];
	//printf("A%d = %f\n",c,diff_seg[2][c]);
}

average_unity/=3;


//printf("%f", diff_seg[0][2]);
/*
printf("Ua:%f;\n", average_unity);
printf("Ud:%f;\n", deviation_unity);
*/


if(deviation_unity != (double)0)
{
	//printf(";Df:%f", average_unity/deviation_unity);
	//printf(";C:%f;\n", (average_unity/deviation_unity));
	//printf(";Oc:%f;", (average_unity/deviation_unity)*diff_seg[0][2]);

	//printf("%f", ((average_unity/deviation_unity)*diff_seg[0][2])*diff_seg[4][0]);
	//printf("%f", ((average_unity/deviation_unity)*diff_seg[0][2]));
	printf("%f", (diff_seg[3][0]+diff_seg[3][1]+diff_seg[3][2])/3);
}
else
{
	//printf(";Df:0");
	//printf(";Oc:%f", (average_unity*deviation_unity)*diff_seg[0][2]);
}



if(debug == 1)
{
	write_bmp_color(image_out, path_fft_x);


//printf("\nGENERAL      DISTANCE: %f\n", diff_general);
//printf("\nSmax: %f\n", diff_seg[0]);
//printf("Smin: %f\n", diff_seg[1]);
//printf("\nSEGMENTATION DISTANCE: %f\n", diff_seg[2]);
//printf("\nET: %fs\n", ((double)clock() - start) / CLOCKS_PER_SEC);

}




	return EXIT_SUCCESS;
}

/*
 * STANDARD FUNCTION
 */

void read_image(unsigned char image_out[2][3][Y_SIZE][X_SIZE], int lim)
{
	char *delim3 = "|";
	char *delim2 = "#";
	char *delim1 = ":";
	char *rtn;
	char *ptn;
	char *btn;
	char *rest;
	char *xest;
	char *yest;
	int xi = 0;
	int yi = 0;
	int ii = 0;
	int ci = 0;
	long fslen = (((lim*lim)*3)*3+ ((lim*lim)*3)*3) *2+lim;
	char *c = malloc(fslen); // allocate memory
	fgets(c, fslen, stdin);

	if((c != NULL)&& !ferror(stdin))
	{
	btn = strtok_r(c, delim3,&yest);
	while(btn != NULL)
	{
		ptn = strtok_r(btn, delim2,&xest);
		while(ptn != NULL)
			{
				rtn = strtok_r(ptn, delim1,&rest);
				while(rtn != NULL)
				{
					if(yi <= lim)
					{
						image_out[ii][ci][xi][yi] = atol(rtn);
						xi++;
					}
					if(xi >= lim)
					{
						yi++;
						xi = 0;
					}
					rtn = strtok_r(NULL, delim1, &rest);
				}
				ptn = strtok_r(NULL, delim2, &xest);
				ci++;
				xi = 0;
				yi = 0;
			}
			btn = strtok_r(NULL, delim3, &yest);
			ii++;
			xi = 0;
			yi = 0;
			ci = 0;
		}
		xi = 0;
		yi = 0;
		ci = 0;
		ii = 0;
	}
	free(c);
}
/*
void sora_decode(ImageData *img, int lim)
{
	char *delim3 = "|";
	char *delim2 = "#";
	char *delim1 = ":";
	char *rtn;
	char *ptn;
	char *btn;
	char *rest;
	char *xest;
	char *yest;
	int xi = 0;
	int yi = 0;
	int ii = 0;
	int ci = 0;
	long fslen = (((lim*lim)*3)*3+ ((lim*lim)*3)*3) *2+lim;
	char *c = malloc(fslen); // allocate memory
	unsigned char image_out[2][3][Y_SIZE][X_SIZE];

	fgets(c, fslen, stdin);

	if((c != NULL)&& !ferror(stdin))
	{
	btn = strtok_r(c, delim3,&yest);
	while(btn != NULL)
	{
		ptn = strtok_r(btn, delim2,&xest);
		while(ptn != NULL)
			{
				rtn = strtok_r(ptn, delim1,&rest);
				while(rtn != NULL)
				{
					if(yi <= lim)
					{
						image_out[ii][ci][xi][yi] = atol(rtn);
						//setPixel(img,xi,yi,Pixel *pix)
						xi++;
					}
					if(xi >= lim)
					{
						yi++;
						xi = 0;
					}
					rtn = strtok_r(NULL, delim1, &rest);
				}
				ptn = strtok_r(NULL, delim2, &xest);
				ci++;
				xi = 0;
				yi = 0;
			}
			btn = strtok_r(NULL, delim3, &yest);
			ii++;
			xi = 0;
			yi = 0;
			ci = 0;
		}
		xi = 0;
		yi = 0;
		ci = 0;
		ii = 0;
	}

	//def
	int i;
	Pixel *pix;

	for(i = 0;i < 2;i++)
	{
		for(c = 0;c < 3;c++)
		{

		    pix->r=val;
		    pix->g=val;
		    pix->b=val;

		}
	}

	free(c);
}
*/

char *join(const char* s1, const char* s2)
{
    char* result = malloc(strlen(s1) + strlen(s2) + 1);

    if (result)
    {
        strcpy(result, s1);
        strcat(result, s2);
    }

    return result;
}

char *tochar(int var)
{
	static char str[1024];
	sprintf(str, "%d", var);
	return str;
}

//sprintf(str, "%d", (time((time_t *) NULL)));

/*
 * STANDARD FUNCTION
 */
void convolute(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int b)
{
	int i, j;
	int c = 1;
	for (i = 1; i < Y_SIZE-1; i++)
	{
		for (j = 1; j < X_SIZE-1; j++)
		{
/*
			image_out[i-1][j-1] = image_in[i-1][j-1] * c;
			image_out[i-1][j] = image_in[i-1][j] * c;
			image_out[i-1][j+1] = image_in[i-1][j+1] * c;
			image_out[i][j-1] = image_in[i][j-1] * c;
			image_out[i][j] = image_in[i][j] * c;
			image_out[i][j+1] = image_in[i][j+1] * c;
			image_out[i+1][j-1]= image_in[i+1][j-1] * c;
			image_out[i+1][j] = image_in[i+1][j] * c;
			image_out[i+1][j] = image_in[i+1][j] * c;
*/
			if ((long)sqrt((double)((j-X_SIZE/2)*(j-X_SIZE/2) + (i-Y_SIZE/2)*(i-Y_SIZE/2))) <= b){
			int centre =
			(image_in[i-1][j-1] * c)
			+(image_in[i-1][j] * c)
			+(image_in[i-1][j+1] * c) //row 1

			+(image_in[i][j-1] * c)
			+(image_in[i][j] * c) //(self)
			+(image_in[i][j+1] * c) //row 2

			+(image_in[i+1][j-1] * c)
			+(image_in[i+1][j] * c)
			+(image_in[i+1][j+1] * c); //row 3

			if(centre >= 255)
				centre = 255;
			if(centre <= 0)
				centre = 0;

			image_out[i][j] = centre;
			}
			else
			{
				image_out[i][j] = image_in[i][j];
			}
		}
	}
}

void assign(unsigned char image_out[Y_SIZE][X_SIZE], unsigned char image_in[Y_SIZE][X_SIZE])
{
	int i, j;
	for (i = 0; i < Y_SIZE; i++)
	{
		for (j = 0; j < X_SIZE; j++)
		{
	image_out[i][j] = image_in[i][j];
		}
	}
}

int write_bmp_mono(unsigned char image[Y_SIZE][X_SIZE], char *filename)
{
	long i, j;
	FILE *fp;
	long file_size, width, height;
	unsigned char *image_buf;
	unsigned char header1[54] = {0x42, 0x4d, 0, 0, 0, 0, 0, 0, 0, 0,
		54, 4, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 8, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0};
	unsigned char header2[1024];
	file_size = (long)X_SIZE * (long)Y_SIZE + 54 + 1024;
	header1[2] = (unsigned char)(file_size & 0x000000ff);
	header1[3] = (unsigned char)((file_size >> 8) & 0x000000ff);
	header1[4] = (unsigned char)((file_size >> 16)  & 0x000000ff);
	header1[5] = (unsigned char)((file_size >> 24)  & 0x000000ff);
	width = X_SIZE;
	header1[18] = (unsigned char)(width & 0x000000ff);
	header1[19] = (unsigned char)((width >> 8) & 0x000000ff);
	header1[20] = (unsigned char)((width >> 16) & 0x000000ff);
	header1[21] = (unsigned char)((width >> 24) & 0x000000ff);
	height = Y_SIZE;
	header1[22] = (unsigned char)(height & 0x000000ff);
	header1[23] = (unsigned char)((height >> 8) & 0x000000ff);
	header1[24] = (unsigned char)((height >> 16) & 0x000000ff);
	header1[25] = (unsigned char)((height >> 24) & 0x000000ff);
	for (i= 0; i < 256; i++) {
		header2[i*4  ] = (unsigned char)i;
		header2[i*4+1] = (unsigned char)i;
		header2[i*4+2] = (unsigned char)i;
		header2[i*4+3] = 0;
	}
	image_buf = (unsigned char *)malloc((size_t)X_SIZE*Y_SIZE);
	if (image_buf == NULL) return -1;
	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			*(image_buf + X_SIZE*i + j) = image[Y_SIZE-i-1][j];
	if ((fp = fopen(filename, "wb")) == NULL) return -1;
	fwrite(header1, sizeof(unsigned char), 54, fp);
	fwrite(header2, sizeof(unsigned char), 1024, fp);
	fwrite(image_buf, sizeof(unsigned char), (size_t)(long)X_SIZE*Y_SIZE, fp);
	fclose(fp);
	free(image_buf);
	return 0;
}


int write_bmp_color(unsigned char image[3][Y_SIZE][X_SIZE], char *filename)
{
	long i, j;
	FILE *fp;
	long file_size, width, height;
	unsigned char *image_buf;
	unsigned char header[54] = {0x42, 0x4d, 0, 0, 0, 0, 0, 0, 0, 0,
		54, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0};
	file_size = (long)X_SIZE * (long)Y_SIZE * 3 + 54;
	header[2] = (unsigned char)(file_size & 0x000000ff);
	header[3] = (unsigned char)((file_size >> 8) & 0x000000ff);
	header[4] = (unsigned char)((file_size >> 16)  & 0x000000ff);
	header[5] = (unsigned char)((file_size >> 24)  & 0x000000ff);
	width = X_SIZE;
	header[18] = (unsigned char)(width & 0x000000ff);
	header[19] = (unsigned char)((width >> 8) & 0x000000ff);
	header[20] = (unsigned char)((width >> 16) & 0x000000ff);
	header[21] = (unsigned char)((width >> 24) & 0x000000ff);
	height = Y_SIZE;
	header[22] = (unsigned char)(height & 0x000000ff);
	header[23] = (unsigned char)((height >> 8) & 0x000000ff);
	header[24] = (unsigned char)((height >> 16) & 0x000000ff);
	header[25] = (unsigned char)((height >> 24) & 0x000000ff);
	image_buf = (unsigned char *)malloc((size_t)X_SIZE*Y_SIZE*3);
	if (image_buf == NULL) return -1;
	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++) {
			*(image_buf + 3*(X_SIZE*i + j)    ) = image[2][Y_SIZE-i-1][j];
			*(image_buf + 3*(X_SIZE*i + j) + 1) = image[1][Y_SIZE-i-1][j];
			*(image_buf + 3*(X_SIZE*i + j) + 2) = image[0][Y_SIZE-i-1][j];
		}
	if ((fp = fopen(filename, "wb")) == NULL) return -1;
	fwrite(header, sizeof(unsigned char), 54, fp);
	fwrite(image_buf, sizeof(unsigned char), (size_t)(long)X_SIZE*Y_SIZE*3, fp);
	fclose(fp);
	free(image_buf);
	return 0;
}


int read_bmp_color(unsigned char image[3][Y_SIZE][X_SIZE], char *filename)
{
long i, j;
FILE *fp;
unsigned char *image_buf;
unsigned char header[54];

image_buf = (unsigned char *)malloc((size_t)X_SIZE*Y_SIZE*3);
if (image_buf == NULL) return -1;
if ((fp = fopen(filename, "rb")) == NULL) return -1;
fread(header, sizeof(unsigned char), 54, fp);
fread(image_buf, sizeof(unsigned char), (size_t)(long)X_SIZE*Y_SIZE*3, fp);
fclose(fp);
for (i = 0; i < Y_SIZE; i++)
	for (j = 0; j < X_SIZE; j++) {
		image[0][Y_SIZE-i-1][j] = *(image_buf + 3 * (X_SIZE*i + j) + 2);
		image[1][Y_SIZE-i-1][j] = *(image_buf + 3 * (X_SIZE*i + j) + 1);
		image[2][Y_SIZE-i-1][j] = *(image_buf + 3 * (X_SIZE*i + j)    );
	}
free(image_buf);
return 0;
}


/*
 * DECL FFT
 * SEE Proto.h
 * DELC FFT END
 */

/*--- fft1 ---------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------------------------------*/
int fft1(double a_rl[], double a_im[], int ex, int inv)
{
	int i, length = 1;
	double *sin_tbl;
	double *cos_tbl;
	double *buf;

	for (i = 0; i < ex; i++) length *= 2;
	sin_tbl = (double *)malloc((size_t)length*sizeof(double));
	cos_tbl = (double *)malloc((size_t)length*sizeof(double));
	buf = (double *)malloc((size_t)length*sizeof(double));
	if ((sin_tbl == NULL) || (cos_tbl == NULL) || (buf == NULL)) return -1;
	cstb(length, inv, sin_tbl, cos_tbl);
	fft1core(a_rl, a_im, length, ex, sin_tbl, cos_tbl, buf);
	free(sin_tbl);
	free(cos_tbl);
	return 0;
}

/*--- fft1core -------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------------------------------*/
void fft1core(double a_rl[], double a_im[], int length, int ex, double sin_tbl[], double cos_tbl[], double buf[])
{
	int i, j, k, w, j1, j2;
	int numb, lenb, timb;
	double xr, xi, yr, yi, nrml;

	if (OPT == 1) {
		for (i = 1; i < length; i+=2) {
			a_rl[i] = -a_rl[i];
			a_im[i] = -a_im[i];
		}
	}
	numb = 1;
	lenb = length;
	for (i = 0; i < ex; i++) {
		lenb /= 2;
		timb = 0;
		for (j = 0; j < numb; j++) {
			w = 0;
			for (k = 0; k < lenb; k++) {
				j1 = timb + k;
				j2 = j1 + lenb;
				xr = a_rl[j1];
				xi = a_im[j1];
				yr = a_rl[j2];
				yi = a_im[j2];
				a_rl[j1] = xr + yr;
				a_im[j1] = xi + yi;
				xr = xr - yr;
				xi = xi - yi;
				a_rl[j2] = xr*cos_tbl[w] - xi*sin_tbl[w];
				a_im[j2] = xr*sin_tbl[w] + xi*cos_tbl[w];
				w += numb;
			}
			timb += (2*lenb);
		}
		numb *= 2;
	}
	birv(a_rl, length, ex, buf);
	birv(a_im, length, ex, buf);
	if (OPT == 1) {
		for (i = 1; i < length; i+=2) {
			a_rl[i] = -a_rl[i];
			a_im[i] = -a_im[i];
		}
	}
	nrml = 1.0 / sqrt((double)length);
	for (i = 0; i < length; i++) {
		a_rl[i] *= nrml;
		a_im[i] *= nrml;
	}
}

/*--- cstb --- -------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------*/
void cstb(int length, int inv, double sin_tbl[], double cos_tbl[])
{
	int i;
	double xx, arg;

	xx = -PI* 2.0 / (double)length;
	if (inv < 0) xx = -xx;
	for (i = 0; i < length; i++) {
		arg = i * xx;
		sin_tbl[i] = sin(arg);
		cos_tbl[i] = cos(arg);
	}
}

/*--- birv ------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------*/
void birv(double a[], int length, int ex, double b[])
{
	int	i, ii, k, bit;

	for (i = 0; i < length; i++) {
		for (k = 0, ii=i, bit=0; ; bit<<=1, ii>>=1) {
			bit = (ii & 1) | bit;
			if (++k == ex) break;
		}
		b[i] = a[bit];
	}
	for (i = 0; i < length; i++) a[i] = b[i];
}






/*--- fft2 ------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------------------*/
int fft2 (double a_rl[Y_SIZE][X_SIZE], double a_im[Y_SIZE][X_SIZE], int inv)
{
	double *b_rl;
	double *b_im;
	double *hsin_tbl;
	double *hcos_tbl;
	double *vsin_tbl;
	double *vcos_tbl;
	double *buf_x;
	double *buf_y;
	int i;

	b_rl = (double *)calloc((size_t)X_SIZE*Y_SIZE, sizeof(double));
	b_im = (double *)calloc((size_t)X_SIZE*Y_SIZE, sizeof(double));
	hsin_tbl = (double *)calloc((size_t)X_SIZE, sizeof(double));
	hcos_tbl = (double *)calloc((size_t)X_SIZE, sizeof(double));
	vsin_tbl = (double *)calloc((size_t)Y_SIZE, sizeof(double));
	vcos_tbl = (double *)calloc((size_t)Y_SIZE, sizeof(double));
	buf_x = (double *)malloc((size_t)X_SIZE*sizeof(double));
	buf_y = (double *)malloc((size_t)Y_SIZE*sizeof(double));
	if ((b_rl == NULL) || (b_im == NULL)
		|| (hsin_tbl == NULL) || (hcos_tbl == NULL)
		|| (vsin_tbl == NULL) || (vcos_tbl == NULL)
		|| (buf_x == NULL) || (buf_y == NULL)) {
		return -1;
	}
	cstb(X_SIZE, inv, hsin_tbl, hcos_tbl);
	cstb(Y_SIZE, inv, vsin_tbl, vcos_tbl);

	for (i = 0; i < Y_SIZE; i++) {
		fft1core(&a_rl[(long)i][0], &a_im[(long)i][0],
					X_SIZE, X_EXP, hsin_tbl, hcos_tbl, buf_x);
	}

	rvmtx1((double (*)[X_SIZE])a_rl, (double (*)[X_SIZE])b_rl, X_SIZE, Y_SIZE);
	rvmtx1((double (*)[X_SIZE])a_im, (double (*)[X_SIZE])b_im, X_SIZE, Y_SIZE);

	for (i = 0; i < X_SIZE; i++) {
		fft1core(&b_rl[(long)Y_SIZE*i], &b_im[(long)Y_SIZE*i],
					Y_SIZE, Y_EXP, vsin_tbl, vcos_tbl, buf_y);
	}

	rvmtx2((double (*)[Y_SIZE])b_rl, (double (*)[Y_SIZE])a_rl, X_SIZE, Y_SIZE);
	rvmtx2((double (*)[Y_SIZE])b_im, (double (*)[Y_SIZE])a_im, X_SIZE, Y_SIZE);
	free(b_rl);
	free(b_im);
	free(hsin_tbl);
	free(hcos_tbl);
	free(vsin_tbl);
	free(vcos_tbl);
	return 0;
}

/*--- rvmtx1 ------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------*/
void rvmtx1(double a[Y_SIZE][X_SIZE], double b[X_SIZE][Y_SIZE],
	int xsize, int ysize)
{
	int i, j;

	for (i = 0; i < ysize; i++)
		for (j = 0; j < xsize; j++)
			b[j][i] = a[i][j];
}

/*--- rvmtx2 ----------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------*/
void rvmtx2(double a[X_SIZE][Y_SIZE], double b[Y_SIZE][X_SIZE],
	int xsize, int ysize)
{
	int i, j;

	for (i = 0; i < ysize; i++)
		for (j = 0; j < xsize; j++)
			b[i][j] = a[j][i];
}






/*--- fftimage -----------------------------------------------------------------------------

--------------------------------------------------------------------------------------------------------------------*/
int fftimage(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE])
{
	double *ar;
	double *ai;
	double norm, max;
	double data;
	long i, j;

	ar = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	ai = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	if ((ar == NULL) || (ai == NULL)) return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ar[X_SIZE*i + j] = (double)image_in[i][j];
			ai[X_SIZE*i + j] = 0.0;
		}
	}

	if (fft2((double (*)[X_SIZE])ar, (double (*)[X_SIZE])ai, 1) == -1) return -1;

	max = 0;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			norm = ar[X_SIZE*i + j]*ar[X_SIZE*i + j]
			     + ai[X_SIZE*i + j]*ai[X_SIZE*i + j];
			if (norm != 0.0) norm = log(norm) / 2.0;
			else norm = 0.0;
			ar[X_SIZE*i + j] = norm;
			if (norm > max) max = norm;
		}
	}
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ar[X_SIZE*i + j] = ar[X_SIZE*i + j]*255 / max;
		}
	}

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			data = ar[X_SIZE*i + j];
			if (data > 255) data = 255;
			if (data <   0) data = 0;
			image_out[i][j] = (unsigned char)data;
		}
	}
	free(ar);
	free(ai);
	return 0;
}




//List 11.4

/*--- fftfilter -----------------------------------------------------

----------------------------------------------------------------------------------------------------*/
int fftfilter(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int a, int b)
{
	double *ar;
	double *ai;
	double *ff;
	double data;
	long i, j, circ;

	ar = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	ai = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	ff = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	if ((ar == NULL) || (ai == NULL) || (ff == NULL)) return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ar[X_SIZE*i + j] = (double)image_in[i][j];
			ai[X_SIZE*i + j] = 0.0;
		}
	}

	if (fft2((double (*)[X_SIZE])ar, (double (*)[X_SIZE])ai, 1) == -1)
		return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for(j = 0; j < X_SIZE; j++) {
			data = (double)((j-X_SIZE/2)*(j-X_SIZE/2)
				+ (i-Y_SIZE/2)*(i-Y_SIZE/2));
			circ = (long)sqrt(data);
			if ((circ >= a) && (circ <= b))
				ff[X_SIZE*i + j] = 1.0;
			else
				ff[X_SIZE*i + j] = 0.0;
		}
	}

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ar[X_SIZE*i + j] *= ff[X_SIZE*i + j];
			ai[X_SIZE*i + j] *= ff[X_SIZE*i + j];
		}
	}

	if (fft2((double (*)[X_SIZE])ar, (double (*)[X_SIZE])ai, -1) == -1)
		return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			data = ar[X_SIZE*i + j];
			if (data > 255) data = 255;
			if (data <   0) data = 0;
			image_out[i][j] = (unsigned char)data;
		}
	}
	free(ar);
	free(ai);
	free(ff);
	return 0;
}
/*
 * FFT FILTER WITH ORIGINAL OUTPUT
 */
int fftfilter_ff(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int a, int b)
{
	double *ar;
	double *ai;
	double *ff;
	double data;
	long i, j, circ;

	ar = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	ai = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	ff = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	if ((ar == NULL) || (ai == NULL) || (ff == NULL)) return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ar[X_SIZE*i + j] = (double)image_in[i][j];
			ai[X_SIZE*i + j] = 0.0;
		}
	}

	if (fft2((double (*)[X_SIZE])ar, (double (*)[X_SIZE])ai, 1) == -1)
		return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for(j = 0; j < X_SIZE; j++) {
			data = (double)((j-X_SIZE/2)*(j-X_SIZE/2)
				+ (i-Y_SIZE/2)*(i-Y_SIZE/2));
			circ = (long)sqrt(data);
			if ((circ >= a) && (circ <= b))
				ff[X_SIZE*i + j] = 1.0;
			else
				ff[X_SIZE*i + j] = 0.0;
		}
	}

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ar[X_SIZE*i + j] *= ff[X_SIZE*i + j];
			ai[X_SIZE*i + j] *= ff[X_SIZE*i + j];
		}
	}
// START NEW SECTION
//THIS SECTION IS TAKEN FROM
//DECL
	double norm, max;
//DECL END

	max = 0;
		for (i = 0; i < Y_SIZE; i++) {
			for (j = 0; j < X_SIZE; j++) {
				norm = ar[X_SIZE*i + j]*ar[X_SIZE*i + j]
				     + ai[X_SIZE*i + j]*ai[X_SIZE*i + j];
				if (norm != 0.0) norm = log(norm) / 2.0;
				else norm = 0.0;
				ar[X_SIZE*i + j] = norm;
				if (norm > max) max = norm;
			}
		}
		for (i = 0; i < Y_SIZE; i++) {
			for (j = 0; j < X_SIZE; j++) {
				ar[X_SIZE*i + j] = ar[X_SIZE*i + j]*255 / max;
			}
		}

		for (i = 0; i < Y_SIZE; i++) {
			for (j = 0; j < X_SIZE; j++) {
				data = ar[X_SIZE*i + j];
				if (data > 255) data = 255;
				if (data <   0) data = 0;
				image_out[i][j] = (unsigned char)data;
			}
		}



// END NEW SECTION





/*
	if (fft2((double (*)[X_SIZE])ar, (double (*)[X_SIZE])ai, -1) == -1)
		return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			data = ar[X_SIZE*i + j];
			if (data > 255) data = 255;
			if (data <   0) data = 0;
			image_out[i][j] = (unsigned char)data;
		}
	}
*/
	free(ar);
	free(ai);
	free(ff);
	return 0;
}

int fftimage_ff(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int a, int b)
{
	double *ar;
	double *ai;
	double *ff;
	double data;
	long i, j, circ;

	ar = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	ai = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	ff = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	if ((ar == NULL) || (ai == NULL) || (ff == NULL)) return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ar[X_SIZE*i + j] = (double)image_in[i][j];
			ai[X_SIZE*i + j] = (double)image_in[i][j];
		}
	}


	if (fft2((double (*)[X_SIZE])ar, (double (*)[X_SIZE])ai, -1) == -1)
		return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			data = ar[X_SIZE*i + j];
			if (data > 255) data = 255;
			if (data <   0) data = 0;
			image_out[i][j] = (unsigned char)data;
		}
	}

	free(ar);
	free(ai);
	free(ff);
	return 0;
}

int fftfeature_ff(unsigned char image_a[Y_SIZE][X_SIZE],
	unsigned char image_b[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int a, int b)
{
	double *ara;
	double *aia;

	double *arb;
	double *aib;

	double *ff;
	double data;
	long i, j, circ;

	ara = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	aia = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));

	arb = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	aib = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));


	ff = (double *)malloc((size_t)Y_SIZE*X_SIZE*sizeof(double));
	if ((ara == NULL) || (aia == NULL) || (ff == NULL)) return -1;
	if ((arb == NULL) || (aib == NULL) || (ff == NULL)) return -1;



	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			ara[X_SIZE*i + j] = (double)image_a[i][j];
			aia[X_SIZE*i + j] = 0.0;

			arb[X_SIZE*i + j] = (double)image_b[i][j];
			aib[X_SIZE*i + j] = 0.0;
		}
	}

	if (fft2((double (*)[X_SIZE])ara, (double (*)[X_SIZE])aia, 1) == -1) return -1;
	if (fft2((double (*)[X_SIZE])arb, (double (*)[X_SIZE])aib, 1) == -1) return -1;

	for (i = 0; i < Y_SIZE; i++) {
		for(j = 0; j < X_SIZE; j++) {
			data = (double)((j-X_SIZE/2)*(j-X_SIZE/2)
				+ (i-Y_SIZE/2)*(i-Y_SIZE/2));
			circ = (long)sqrt(data);
			if ((circ >= a) && (circ <= b))
				ff[X_SIZE*i + j] = 1.0;
			else
				ff[X_SIZE*i + j] = 0.0;
		}
	}



	double norm_a, norm_b, max_a, max_b;

	for (i = 0; i < Y_SIZE; i++) {
					for (j = 0; j < X_SIZE; j++) {
						norm_a = ara[X_SIZE*i + j]*ara[X_SIZE*i + j]
						     + aia[X_SIZE*i + j]*aia[X_SIZE*i + j];

						norm_b = arb[X_SIZE*i + j]*arb[X_SIZE*i + j]
							 + aib[X_SIZE*i + j]*aib[X_SIZE*i + j];


						if (norm_a != 0.0) norm_a = log(norm_a) / 2.0;
						else norm_a = 0.0;

						if (norm_b != 0.0) norm_b = log(norm_b) / 2.0;
						else norm_b = 0.0;

						//ara[X_SIZE*i + j] = norm_a;
						//arb[X_SIZE*i + j] = norm_b;


						if (norm_a > max_a) max_a = norm_a;
						if (norm_b > max_b) max_b = norm_b;
					}
				}

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			/*
			ara[X_SIZE*i + j] *= ff[X_SIZE*i + j];
			aia[X_SIZE*i + j] *= ff[X_SIZE*i + j];

			arb[X_SIZE*i + j] *= ff[X_SIZE*i + j];
			aib[X_SIZE*i + j] *= ff[X_SIZE*i + j];
			*/

			ara[X_SIZE*i + j] = ara[X_SIZE*i + j] * arb[X_SIZE*i + j] /norm_b;
			aia[X_SIZE*i + j] = aia[X_SIZE*i + j] * aib[X_SIZE*i + j] /norm_b;
	}
	}
	/*
	// START NEW SECTION
	//THIS SECTION IS TAKEN FROM
	//DECL

	//DECL END

		max = 0;
			for (i = 0; i < Y_SIZE; i++) {
				for (j = 0; j < X_SIZE; j++) {
					norm = ara[X_SIZE*i + j]*ara[X_SIZE*i + j]
					     + aia[X_SIZE*i + j]*aia[X_SIZE*i + j];
					if (norm != 0.0) norm = log(norm) / 2.0;
					else norm = 0.0;
					ara[X_SIZE*i + j] = norm;
					if (norm > max) max = norm;
				}
			}
			for (i = 0; i < Y_SIZE; i++) {
				for (j = 0; j < X_SIZE; j++) {
					ara[X_SIZE*i + j] = ara[X_SIZE*i + j]*255 / max;
				}
			}

			for (i = 0; i < Y_SIZE; i++) {
				for (j = 0; j < X_SIZE; j++) {
					data = ara[X_SIZE*i + j];
					if (data > 255) data = 255;
					if (data <   0) data = 0;
					image_out[i][j] = (unsigned char)data;
				}
			}



	// END NEW SECTION



	 */


		if (fft2((double (*)[X_SIZE])ara, (double (*)[X_SIZE])aia, -1) == -1)
			return -1;

		for (i = 0; i < Y_SIZE; i++) {
			for (j = 0; j < X_SIZE; j++) {
				data = ara[X_SIZE*i + j];
				if (data > 255) data = 255;
				if (data <   0) data = 0;
				image_out[i][j] = (unsigned char)data;
			}
		}

		free(ara);
		free(aia);
		free(arb);
		free(aib);
		free(ff);
		return 0;
}

/*
 * RGB TO YSH / YSH TO RGB
 */
void rgb_to_ysh(unsigned char image_in_rgb[3][Y_SIZE][X_SIZE],
	int image_out_ysh[3][Y_SIZE][X_SIZE])
{
	int    i, j;
	double r, g, b, y, cb, cr, s, h;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			r = (double)image_in_rgb[0][i][j];
			g = (double)image_in_rgb[1][i][j];
			b = (double)image_in_rgb[2][i][j];
			y  = 0.2126 * r + 0.7152 * g + 0.0722 * b;
			cb = (b - y) / 1.8556;
			cr = (r - y) / 1.5748;
			s  = sqrt(cb * cb + cr * cr);
			if (s != 0) h  = atan2(cr, cb) * 180.0 / PI;
			else        h = 0;
			image_out_ysh[0][i][j] = (int)y;
			image_out_ysh[1][i][j] = (int)s;
			image_out_ysh[2][i][j] = (int)h;
		}
	}
}


void ysh_to_rgb(int image_in_ysh[3][Y_SIZE][X_SIZE],
	unsigned char image_out_rgb[3][Y_SIZE][X_SIZE])
{
	int    i, j;
	double r, g, b, y, cb, cr, rad;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			y  = image_in_ysh[0][i][j];
			rad = (double)(PI * image_in_ysh[2][i][j] / 180.0);
			cb = image_in_ysh[1][i][j] * cos(rad);
			cr = image_in_ysh[1][i][j] * sin(rad);
			r = y + 1.5748 * cr;
			b = y + 1.8556 * cb;
			g = (y - 0.2126 * r - 0.0722 * b) / 0.7152;
			if (r <   0) r =   0;
			if (r > 255) r = 255;
			if (g <   0) g =   0;
			if (g > 255) g = 255;
			if (b <   0) b =   0;
			if (b > 255) b = 255;
			image_out_rgb[0][i][j] = (unsigned char)r;
			image_out_rgb[1][i][j] = (unsigned char)g;
			image_out_rgb[2][i][j] = (unsigned char)b;
		}
	}
}
/*
 * RGB TO YSH / YSH TO RGB END
 */
/*
 * THRESHOLD
 */
// get brightness
void y_image(int image_in_y[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE])
{
	int i, j, d;

	for (i = 0; i < Y_SIZE; i++){
		for (j = 0; j < X_SIZE; j++ ){
			d = image_in_y[i][j];
			if (d <   0) d =   0;
			if (d > 255) d = 255;
			image_out[i][j] = (unsigned char)d;
		}
	}
}

int threshdiscrim(long hist[256], double disparity)
{
	int i, k;
	double n0, n1, n2, m0, m1, m2;
	double v[256], vmax, v0;

	n0 = 0.0;
	m0 = 0.0;
	for (i = 0; i < 256; i++) {
		n0 += hist[i];
		m0 += i * hist[i];
	}
	if (n0 == 0.0) m0 = 0.0;
	else m0 /= n0;
	v0 = 0.0;
	for (i = 0; i < 256; i++) v0 += hist[i] * (i - m0) * (i - m0) / n0;
	for (k = 0; k < 256; k++) {
		n1 = 0.0;
		m1 = 0.0;
		for (i = 0; i < k; i++) {
			n1 += hist[i];
			m1 += i * hist[i];
		}
		if (n1 == 0.0) m1 = 0.0;
		else m1 /= n1;
		n2 = 0.0;
		m2 = 0.0;
		for (i = k; i < 256; i++) {
			n2 += hist[i];
			m2 += i * hist[i];
		}
		if (n2 == 0.0) m2 = 0.0;
		else m2 /= n2;
		v[k] = (n1 * (m1 - m0) * (m1 - m0) + n2 * (m2 - m0) * (m2 - m0)) / n0;
	}
	vmax = 0.0;
	for (i = 0; i < 256; i++) {
		if (vmax <= v[i]) {
			vmax = v[i];
			k = i;
		}
	}
	if (v0 == 0) return 0;
	if ((vmax / v0) >= disparity) return k;
	else return 0;
}




/*--- threshold_dynamic ----------------------------------------------------------------------------------*/
void threshold_dynamic(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int type)
{
	int	 i, j, k, m, n, m1, m2, n1, n2, s, t;
	int  thm[DIV+1][DIV+1];
	long hist[256];
	int  thresh;
	double p, q;

	for (i = 0; i <= DIV; i++) {
		for (j = 0; j <= DIV; j++) {
			thm[i][j] = 0;
		}
	}

	for (i = 0; i <= DIV; i++) {
		for (j = 0; j <= DIV; j++) {
			for (k = 0; k < 256; k++) hist[k] = 0;
			if (i != 0) m1 = -YS;
			else m1 = 0;
			if (i != DIV) m2 = YS;
			else m2 = 0;
			if (j != 0) n1 = -XS;
			else n1 = 0;
			if (j != DIV) n2 = XS;
			else n2 = 0;
			for (m = m1; m < m2; m++) {
				for (n = n1; n < n2; n++) {
					k = image_in[i*YS+m][j*XS+n];
					hist[k]++;
				}
			}
			thm[i][j] = threshdiscrim(hist, DTH);
		}
	}

	for (i = 0; i <= DIV; i++) {
		for (j = 0; j <= DIV; j++) {
			if (thm[i][j] <= 0) {
				for (k = 1; k < DIV; k++) {
					s = 0;
					t = 0;
					m1 = i - k;
					m2 = i + k;
					n1 = j - k;
					n2 = j + k;
					if (m1 <   0) m1 = 0;
					if (m2 > DIV) m2 = DIV;
					if (n1 <   0) n1 = 0;
					if (n2 > DIV) n2 = DIV;
					for (m = m1; m <= m2; m++) {
						for (n = n1; n <= n2; n++) {
							if (thm[m][n] > 0) {
								s += 1 / k;
							    t += thm[m][n] / k;
							}
						}
					}
					if (s >= 4) {
						thm[i][j] = t / s;
						break;
					}
				}
			}
		}
	}

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			m = i / YS;
			n = j / XS;
			q = (double)(i % YS) / YS;
			p = (double)(j % XS) / XS;
			thresh = (int)((1.0-q)*((1.0-p)*thm[m  ][n  ]
			                            + p*thm[m  ][n+1])
			                   + q*((1.0-p)*thm[m+1][n  ]
			                            + p*thm[m+1][n+1]));
			switch (type){
				case 2:
					if ((int)image_in[i][j] <= thresh)
						image_out[i][j] = HIGH;
					else
						image_out[i][j] =  LOW;
					break;
				default:
					if ((int)image_in[i][j] >= thresh)
						image_out[i][j] = HIGH;
					else
						image_out[i][j] =  LOW;
					break;
			}
		}
	}
}
void threshold_discrim(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int type)
{
	int	i, j;
	int thresh;
	long hist[256];

	histgram(image_in, hist);
	thresh = threshdiscrim(hist, 0.0);
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			switch (type){
				case 2:
					if ((int)image_in[i][j] <= thresh) image_out[i][j] = HIGH;
					else                               image_out[i][j] =  LOW;
					break;
				default:
					if ((int)image_in[i][j] >= thresh) image_out[i][j] = HIGH;
					else                               image_out[i][j] =  LOW;
					break;
			}
		}
	}
}
void histgram(unsigned char image_in[Y_SIZE][X_SIZE], long hist[256])
{
	int	i, j, n;

	for (n = 0; n < 256; n++) hist[n] = 0;
	for (i = 0; i < Y_SIZE; i++) {
		for ( j = 0; j < X_SIZE; j++) {
			n = image_in[i][j];
			hist[n]++;
		}
	}
}

/*
 * THRESHOLD END
 */
/*
 * FEATURE
 */





/*--- labeling --------------------*/
int labeling(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_label[Y_SIZE][X_SIZE])
{
	int	i, j, label;

	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			image_label[i][j] = image_in[i][j];
	label = L_BASE;
	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++) {
			if (image_label[i][j] == HIGH) {
				if (label >= HIGH) return -1;
				labelset(image_label, j, i, label); label++;
			}
	    }
	return label - L_BASE;
}

/*--- labelset -----------*/
void labelset(unsigned char	image[Y_SIZE][X_SIZE], int xs, int ys, int label)
{
	int	i, j, cnt, im, ip, jm, jp;

	image[ys][xs] = label;
	for (;;) {
		cnt = 0;
		for (i = 0; i < Y_SIZE; i++)
			for (j = 0; j < X_SIZE; j++)
				if (image[i][j] == label) {
					im = i-1; ip = i+1; jm = j-1; jp = j+1;
					if (im < 0) im = 0; if (ip >= Y_SIZE) ip = Y_SIZE-1;
					if (jm < 0) jm = 0; if (jp >= X_SIZE) jp = X_SIZE-1;
					if (image[i ][jp] == HIGH) {
						image[i ][jp] = label; cnt++;
					}
					if (image[im][jp] == HIGH) {
						image[im][jp] = label; cnt++;
					}
					if (image[im][j ] == HIGH) {
						image[im][j ] = label; cnt++;
					}
					if (image[im][jm] == HIGH) {
						image[im][jm] = label; cnt++;
					}
					if (image[i ][jm] == HIGH) {
						image[i ][jm] = label; cnt++;
					}
					if (image[ip][jm] == HIGH) {
						image[ip][jm] = label; cnt++;
					}
					if (image[ip][j ] == HIGH) {
						image[ip][j ] = label; cnt++;
					}
					if (image[ip][jp] == HIGH) {
						image[ip][jp] = label; cnt++;
					}
				}
		if (cnt == 0) break;
	}
}





/*--- features -------------------------------*/
void features(unsigned char	image_label_in[Y_SIZE][X_SIZE],
	unsigned char image_label_out[Y_SIZE][X_SIZE],
	int cnt, double size[], double length[], double ratio[], char *buf)
{
	int		i, j, center_x, center_y;
	double	l;
	int		posi, m;

	posi = 0;
	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			image_label_out[i][j] = image_label_in[i][j];
	//m = sprintf(buf, " no     area      circum       round     grav(x,y)\n");
	m = sprintf(buf, " no;area;circum;round;grav(x,y)|\n");
	posi += m;
	for (i = 0; i < cnt; i++) {
		size[i] = calc_size(image_label_out, i+L_BASE,
			&center_x, &center_y);
		l = calc_length(image_label_out, i+L_BASE);
		ratio[i] = 4*PI*size[i]/(l*l);
		image_label_out[center_y][center_x] = HIGH;
		//m = sprintf(&buf[posi], "%3d   %6d   %8.2f   %8.4f     (%3d,%3d)\n",
		//	i, (int)size[i], l, ratio[i], center_x, center_y);
		m = sprintf(&buf[posi], "%3d;%6d;%8.2f;%8.4f;(%3d,%3d)|\n",
					i, (int)size[i], l, ratio[i], center_x, center_y);
		posi += m;
	}
}
void features_compare(unsigned char	image_label_in_a[Y_SIZE][X_SIZE],
	unsigned char	image_label_in_b[Y_SIZE][X_SIZE],
	unsigned char image_label_out_a[Y_SIZE][X_SIZE],
	unsigned char image_label_out_b[Y_SIZE][X_SIZE],
	int cnt_a, int cnt_b, double size[], double length[],
	double ratio[], char *buf, double *rslt, int label[], int *cc, double threshold_ff)
{
	unsigned char image_and[Y_SIZE][X_SIZE];
	unsigned char image_buf_b[Y_SIZE][X_SIZE];
	unsigned char image_buf_a[Y_SIZE][X_SIZE];
	int	i, j;
	int centre_ax[cnt_a], centre_ay[cnt_a];
	int centre_bx[cnt_b], centre_by[cnt_b];
	int a, b;
	//int am, aa;
	//int bm, ba;
	//int cnt_array;
	double d[cnt_a][cnt_b];
	double dr[cnt_a][cnt_b];
	double dl[cnt_a][cnt_b];
	/*
	double stddev[cnt_a];
	double stddev_r[cnt_a];
	double stddev_l[cnt_a];

	double avg[cnt_a];
	double avg_r[cnt_a];
	double avg_l[cnt_a];
	*/

	double cnt_data = 0;
	int cnt_label = 0;
	int label_a[cnt_a];
	int label_b[cnt_b];
	//float total_d;
	int size_a[HIGH];
	int size_b[HIGH];
	double ratio_a[HIGH];
	double ratio_b[HIGH];

	int length_a[HIGH];
	int length_b[HIGH];
	*cc = 0;
	//double thresh = 0.4;
	//double	l[cnt_a][cnt_b];
	int	posi, m;

	posi = 0;
	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			image_label_out_a[i][j] = image_label_in_a[i][j];
	for (i = 0; i < Y_SIZE; i++)
			for (j = 0; j < Y_SIZE; j++)
				image_label_out_b[i][j] = image_label_in_b[i][j];

	m =sprintf(buf, "START=>\n");
	posi += m;

	//DIRECTION : A SEARCH IN B; sample*A => matrix[B]
	for(a = 0; a < cnt_a; a++)
	{
		label_a[a] = a;
		size_a[a] = calc_size(image_label_in_a, a+L_BASE, &centre_ax[a], &centre_ay[a]);
		length_a[a] = calc_length(image_label_in_a, a+L_BASE);
		ratio_a[a] = 4*PI*size_a[a]/(length_a[a]*length_a[a]);

		m=sprintf(&buf[posi], "SECTION A%d=>\n", a);posi += m;
		m=sprintf(&buf[posi], "%6d;|\n", size_a[a]);posi += m;
		for(b = 0; b < cnt_b; b++)
		{
			label_b[b] = b;
			//m=sprintf(&buf[posi], "SECTION B%d=>\n", b);posi += m;
			size_b[b] = calc_size(image_label_in_b, b+L_BASE, &centre_bx[b], &centre_by[b]);
			length_b[b] = calc_length(image_label_in_b, b+L_BASE);
			ratio_b[b] = 4*PI*size_b[b]/(length_b[b]*length_b[b]);

			//RATIO DISTANCE
			if(ratio_a[a] > ratio_b[b])
				{dr[a][b] = (double)ratio_b[b]/(double)ratio_a[a];}
			else if(ratio_a[a] < ratio_b[b])
				{dr[a][b] = (double)ratio_a[a]/(double)ratio_b[b];}
			else if(ratio_a[a] == ratio_b[b])
				{dr[a][b] = 1;}

			//LENGTH DISTANCE
			if(length_a[a] > length_b[b])
				{dl[a][b] = (double)length_b[b]/(double)length_a[a];}
			else if(length_a[a] < length_b[b])
				{dl[a][b] = (double)length_a[a]/(double)length_b[b];}
			else if(length_a[a] == length_b[b])
				{dl[a][b] = 1;}
			//printf("LENGTH A %d =>%f\n", a,dl[a][b]);
			//printf("RATIO A %d =>%f\n", a,dr[a][b]);

			//SHAPE DISTANCE
			d[a][b] = calc_distance(
					image_label_in_a, image_label_in_b, image_and,
					image_buf_a,
					image_buf_b,
					size_a[a], size_b[b],
				centre_ax[a], centre_ay[a], centre_bx[b], centre_by[b], label_a[a], label_b[b]);
			//printf("\nCURRENT DISTANCE VALUE=> %f FROM %dA, %dB\n", d[a][b], a,b);
			//l[a][b] = abs(size_b[b]-size_a[a]);
			m=sprintf(&buf[posi], ">DISTANCE:%d; => %f;<\n",b, d[a][b]);posi += m;
		}

		//END Ax[DATA] => Ax $$ Bx
		//SORT Ax[DATA] array

		qsort(d[a], cnt_b, sizeof(double), cmp);
		qsort(dr[a], cnt_b, sizeof(double), cmp);
		qsort(dl[a], cnt_b, sizeof(double), cmp);
		/*
		int z;
		printf("CURRENT ARRAY VALUE=> {");
		for(z = 0; z < cnt_b; z++)
			printf("%f, ", d[a][z]);
		printf("}\n");
		printf("CURRENT ARRAY VALUE=> {");
		for(z = 0; z < cnt_b; z++)
			printf("%f, ", dl[a][z]);
		printf("}\n");
		printf("CURRENT ARRAY VALUE=> {");
		for(z = 0; z < cnt_b; z++)
			printf("%f, ", dr[a][z]);
		printf("}\n");
		*/
		//stddev[a] = stat_stddev(d[a], cnt_b);
		//stddev_r[a] = stat_stddev(dr[a], cnt_b);
		//stddev_l[a] = stat_stddev(dl[a], cnt_b);
		//avg[a] = stat_avg(d[a], cnt_b);
		//avg_r[a] = stat_avg(dr[a], cnt_b);
		//avg_l[a] = stat_avg(dl[a], cnt_b);
		//printf("STDDEV A&B %d =>%f\n", a,stddev_l[a]);
		//qsort(&l[a][b], cnt_b, sizeof(double), sort);
					//for(bm = 0; bm < cnt_b; bm++)
						//ba += d[a][bm];
					//ba /= cnt_b;//CALC AVERAGE
					//PEAK DETECT/UNIQUE DETECT
					//CHECK COMPARISON RESULT
					m=sprintf(&buf[posi], ">VALUE => H=%f, L=%f\n",d[a][0], d[a][cnt_b-1]);posi += m;
					m=sprintf(&buf[posi], ">VALUE => H-L =  %f\n",d[a][0]- d[a][cnt_b-1]);posi += m;

					//printf("STDDEV A&B %d =>%f AVG => %f HIGH => %f LOW => %f\n", a, stddev[a], avg[a], d[a][0],d[a][cnt_b-1]);
					/*
					if(
							((d[a][0]) >= threshold_ff)

							||

							//((dl[a][0]-dl[a][0+1]) <= threshold_ff)

							//||

							((dr[a][0]-dr[a][0+1]) <= threshold_ff)

					)
					{
					*/
						//printf("PEAK!\n");


						m=sprintf(&buf[posi], ">PEAK DATA LOGIC => %f\n", (d[a][0]-d[a][0+1])); posi += m;
						m=sprintf(&buf[posi], ">PEAK DATA LENGTH=> %f\n", (dl[a][0]-dl[a][0+1])); posi += m;
						m=sprintf(&buf[posi], ">PEAK DATA RATIO => %f\n", (dr[a][0]-dr[a][0+1])); posi += m;
						double d_dif[cnt_b];
						double dl_dif[cnt_b];
						double dr_dif[cnt_b];
						//printf("#########STAT DATA %d###PREPARE\n", a);
						if((cnt_b-1) > 2)//check value COUNT for stat curve analysis,
						{
							for(i = 0; i < cnt_b-1; i++)
							{
								d_dif[i] = fabs(d[a][i]-d[a][i+1]);
								dl_dif[i] = fabs(dl[a][i]-dl[a][i+1]);
								dr_dif[i] = fabs(dr[a][i]-dr[a][i+1]);
							}
							//double d_stddev = stat_stddev(d_dif, cnt_b-1);
							//double dl_stddev = stat_stddev(dl_dif, cnt_b-1);
							//double dr_stddev = stat_stddev(dr_dif, cnt_b-1);


							double imagine_point = (double)(d[a][0] + d[a][cnt_b]) / 2;
							double real_point = d[a][((int)round((double)cnt_b/(double)2))];
							//1-value -> opposite

							double stat_variation;
							if(imagine_point > real_point)
							{
								stat_variation = (double)fabs((double)imagine_point - (double)real_point)/(double)fmax(real_point, imagine_point);
							}
							else if(imagine_point < real_point)
							{
								stat_variation = 1-fabs(imagine_point - real_point)/fmax(real_point, imagine_point);
							}
							else
							{
								stat_variation = 0;
							}
							//second way, not recommended:
							//double stat_variation = 1 - fabs(fmin(imagine_point / real_point, real_point / imagine_point));
							//
							double stat_ana = pow(stat_variation, pow((double)10*((double)1-stat_variation),pow((double)1/((double)1-stat_variation), (double)1/((double)1-stat_variation))));
							//y=x^{10|_cdot_({1-x})^{|_frac_{{1};{({1-x})}}^{|_frac_{{1};{({1-x})}}}}} <- detection equation
							/*
							printf("#########STAT DATA %d#####START\n", a);
							printf(">POINT 0       => %f\n", d[a][0]);
							printf(">POINT n       => %f\n", d[a][cnt_b]);
							printf(">INDEX 0       => %d\n", 0);
							printf(">INDEX n       => %d\n", cnt_b);
							printf(">INDEX (n+0)/2 => %d\n", ((int)round((double)cnt_b/(double)2)));
							printf(">IMAGINE POINT => %f\n", imagine_point);
							printf(">REAL    POINT => %f\n", real_point);
							printf(">REAL    INDEX => %f\n", real_point);
							// T OR F
							printf(">VARIATION(DIF)=> %f\n", stat_variation);
							printf(">PEAK VALUE    => %f\n", stat_ana);
							printf("#########STAT DATA %d#######END\n", a);
							*/
							cnt_data = cnt_data+stat_ana;
							//GET PEAK DATA; MARK AS FOUND WITH ANA STATS
							cnt_label++;
							//LABEL PEAK SEGMENT
							label[cnt_label] = (a+L_BASE);
							//LABEL CONSTANTS SET


							/*
							double stat_final =
									((double)stat_variation *d_stddev)
									//((double)stat_variation *dl_stddev)
									//((double)stat_variation *dr_stddev)
									;//double ends here...
							*/



							/*
							if(
									imagine_point > real_point
							)
							{
								//VALUE DOES NOT VARY
								//printf("\nSTDDEV LOGIC STAT [CURVE]=> :%f \n", d_stddev);
								//printf("\nSTDDEV LOGIC STAT [CURVE]=> :%f \n", dl_stddev);
								//printf("\nSTDDEV LOGIC STAT [CURVE]=> :%f \n", dr_stddev);
							}
							else if((d[a][0] + d[a][cnt_b])/2 < d[a][((int)round(cnt_b/2))])
							{
								d_stddev = 1-d_stddev;
								dl_stddev = 1-dl_stddev;
								dr_stddev = 1-dr_stddev;
								//printf("\nSTDDEV LOGIC STAT [!CURVE]=> :%f \n", d_stddev);
								//printf("\nSTDDEV LOGIC STAT [!CURVE]=> :%f \n", dl_stddev);
								//printf("\nSTDDEV LOGIC STAT [!CURVE]=> :%f \n", dr_stddev);
							}
							*/

							//cnt_data = cnt_data + (d_stddev+dl_stddev+dr_stddev);

							//printf("MAKING LABEL PROG  %d=> %d\n",cnt_data , label[cnt_data]);
							m=sprintf(&buf[posi], ">PEAK SECTION => %d\n", a);posi += m;
							//CHECK IF ONE LABEL THRESHOLD HAS LARGER VALUE THAN OTHERS IN RESULT(PEAK TARGET DATA):3
							//IF UNDETECTED(GARBAGE DATA), MAIN CNT DECREMENT; BEFORE END, CHECK CNT VALUE,
							//IF CNT = ORIGINAL CNT; THEN => MATRIX A = MATRIX B
							//IF CNT < ORIGINAL CNT *BY VALUE*; THEN => MATRIX A < MATRIX B BY *VALUE*
							//IF CNT = 0; THEN => MATRIX A != MATRIX B
						}
						//printf("#########STAT DATA %d##COMPLETE\n\n", a);
					//}
					/*else
					{

					   //printf("PEAK!\n");
						cnt_data++;//GET PEAK DATA; MARK AS FOUND....
						//TODO FIX CNT_B = 0 => NaN or inf BUG!
						cnt_label++;//LABEL PEAK SEGMENT
						label[cnt_label] = (a+L_BASE);//LABEL
						m=sprintf(&buf[posi], ">PEAK DATA LOGIC => %f\n", (d[a][0]-d[a][0+1]));posi += m;
						m=sprintf(&buf[posi], ">PEAK DATA LENGTH=> %f\n", (dl[a][0]-dl[a][0+1]));posi += m;
						m=sprintf(&buf[posi], ">PEAK DATA RATIO => %f\n", (dr[a][0]-dr[a][0+1]));posi += m;
						double d_dif[cnt_b];
						double dl_dif[cnt_b];
						double dr_dif[cnt_b];
						//double tmp = d[a][0]-d[a][cnt_b];
						//printf("!%f!", tmp);
						//printf("CNT COUNT B = %d, A = %d\n", cnt_b, cnt_a);
						if((cnt_b-1) > 1)
						{
							for(i = 0; i < cnt_b-1; i++)
							{
								d_dif[i] = d[a][i]-d[a][i+1];
								dl_dif[i] = dl[a][i]-dl[a][i+1];
								dr_dif[i] = dr[a][i]-dr[a][i+1];
							}
							//DEBUG

							//int z;
							//printf("CURRENT STDDEV ARRAY VALUE=> {");
							//for(z = 0; z < cnt_b-1; z++)
							//	printf("%f, ", d_dif[z]);
							//printf("}\n");

							//DEBUG
							double d_stddev = stat_stddev(d_dif, cnt_b-1);
							double dl_stddev = stat_stddev(dl_dif, cnt_b-1);
							double dr_stddev = stat_stddev(dr_dif, cnt_b-1);
							if((d[a][0] + d[a][cnt_b])/2 > d[a][((int)round(cnt_b/2))])
							{
								//VALUE DOES NOT VARY
								//printf("\nSTDDEV LOGIC STAT [CURVE]=> :%f \n", d_stddev);
								//printf("\nSTDDEV LOGIC STAT [CURVE]=> :%f \n", dl_stddev);
								//printf("\nSTDDEV LOGIC STAT [CURVE]=> :%f \n", dr_stddev);
							}
							else if((d[a][0] + d[a][cnt_b])/2 < d[a][((int)round(cnt_b/2))])
							{
								d_stddev = 1-d_stddev;
								dl_stddev = 1-dl_stddev;
								dr_stddev = 1-dr_stddev;
								//printf("\nSTDDEV LOGIC STAT [!CURVE]=> :%f \n", d_stddev);
								//printf("\nSTDDEV LOGIC STAT [!CURVE]=> :%f \n", dl_stddev);
								//printf("\nSTDDEV LOGIC STAT [!CURVE]=> :%f \n", dr_stddev);
							}

							cnt_data =  (d_stddev+dl_stddev+dr_stddev);
							//printf("MAKING LABEL PROG  %d=> %d\n",cnt_data , label[cnt_data]);
							m=sprintf(&buf[posi], ">PEAK SECTION => %d\n", a);posi += m;
							//CHECK IF ONE LABEL THRESHOLD HAS LARGER VALUE THAN OTHERS IN RESULT(PEAK TARGET DATA):3
							//IF UNDETECTED(GARBAGE DATA), MAIN CNT DECREMENT; BEFORE END, CHECK CNT VALUE,
							//IF CNT = ORIGINAL CNT; THEN => MATRIX A = MATRIX B
							//IF CNT < ORIGINAL CNT *BY VALUE*; THEN => MATRIX A < MATRIX B BY *VALUE*
							//IF CNT = 0; THEN => MATRIX A != MATRIX B
						}
					}
					*/
					//HERE, LOOP THROUGH THR CNT MATRIX; SORT THRESHOLD DATA
					//CALCULATE DIFFERENCE BETWEEN ARRAY DAT A AND FIND UNIQUE VALUE(PEAK TARGET DATA):3
					//ALG => FIND UNIQUE DATA(ONE)
					//IF CNT
					//(CNT A)
					//image_label_out_a[center_y][center_x] = HIGH;
	}

	*rslt = (double)cnt_data/(double)cnt_a;

	//int c;
	//for(c = 1; c < cnt_data+1; c++)


	*cc = cnt_data;
	m=sprintf(&buf[posi], "<= CC END HERE RSLT => %f,%f, %d \n", *rslt,cnt_data,cnt_a);
	posi += m;
	// rslt = (double)cnt_data/(double)cnt_a;
	// RETURN cnt_data[cnt_a] ?? REMOVE VOID

}
/*--- calc_size ----*/
double calc_size(unsigned char image_label[Y_SIZE][X_SIZE],
	int label, int *cx, int *cy)
{
	int		i, j;
	double	tx, ty, total;
//
//
	tx = 0; ty = 0; total = 0;
	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			if (image_label[i][j] == label) {
				tx += j; ty += i; total++;
			}

	if (total == 0.0) return 0.0;
	*cx = (int)(tx/total); *cy = (int)(ty/total);
	return total;
}

int cmp(const void *x, const void *y)
{
	double xx = *(double*)x, yy = *(double*)y;
	if (xx > yy) return -1;
	if (xx < yy) return 1;
	return 0;
}

int cmp_struct(const void *x, const void *y)
{
    t_array_index *struct_x = (t_array_index *) x;
    t_array_index *struct_y = (t_array_index *) y;
	if (struct_x->value > struct_y->value) return -1;
	if (struct_x->value < struct_y->value) return 1;
	return 0;
}


double stat_stddev(double stat[], int stat_count)
{
	int	i;
	double avg_c = (double)0;
	double stat_c;
	double stddev;
	for(i = 0; i < stat_count; i++)
		avg_c = avg_c + stat[i];
	avg_c = avg_c/stat_count;
	for(i = 0; i < stat_count-1; i++)
		stat_c = pow(stat[i]-avg_c, 2);
	stddev = sqrt((stat_c/stat_count));
	return stddev;
}
double stat_avg(double stat[], int stat_count)
{
	int	i;
	double avg_c = 0;
	for(i = 0; i < stat_count; i++)
		avg_c = avg_c + stat[i];
	if(stat_count != 0)
		avg_c /= stat_count;
	else
		avg_c = (double)0;
	return avg_c;
}
double euclidean_dist(double array[][2], int count)
{
	int i;
	double unity = 0;
	for(i = 0; i < count; i++)
		unity += pow((array[i][0]-array[i][1]), 2);
	return sqrt(unity);
}

double calc_distance(unsigned char image_label_a[Y_SIZE][X_SIZE],
	unsigned char image_label_b[Y_SIZE][X_SIZE],
	unsigned char image_label_and[Y_SIZE][X_SIZE],
	unsigned char image_buf_a[Y_SIZE][X_SIZE],
	unsigned char image_buf_b[Y_SIZE][X_SIZE],
	int size_a, int size_b,
	int centre_ax, int centre_ay,
	int centre_bx, int centre_by,
	int label_a, int label_b)
{
	int total = 0;
	//SHIFT MATRIX TO CENTRE
	int shift_ax, shift_ay;
	int shift_bx, shift_by;
	//SHIFT A TO B
	shift_ax = centre_bx - centre_ax;
	shift_ay = centre_by - centre_ay;
	//SHIFT B TO A
	shift_bx = centre_ax - centre_bx;
	shift_by = centre_ay - centre_by;
	//printf("\n=>SHIFT XY => !");
	//printf("LOGIC INTERSECTION: %d => A(%d, %d); SHIFT => X:%d, Y:%d\n", label_a, centre_ax, centre_ay, shift_ax, shift_ay);
	//printf("LOGIC INTERSECTION: %d => B(%d, %d); SHIFT => X:%d, Y:%d\n", label_b, centre_bx, centre_by, shift_bx, shift_by);

	int	i, j;
	//CENTRE POINT SHIFT END
	//LOOP MATRIX CONTENT

	for (i = 0; i < Y_SIZE; i++)
	{
		for (j = 0; j < X_SIZE; j++)
		{

			if(
					//(shift_bx > shift_ax) && (shift_by > shift_ay)
					1==1
			)
			{
				if(image_label_a[i][j] == (label_a+L_BASE))
				{
					if(
							((i+shift_bx > 0)&&(j+shift_by > 0))        &&
							((i+shift_bx < X_SIZE)&&(j+shift_by < Y_SIZE))
					)
					{
					if(image_label_b[i+shift_bx][j+shift_by] == (label_b+L_BASE))
						{
							total++;

							if((i+shift_bx < 0)||(j+shift_by < 0))
							{
								printf("SHIFT ERROR!:BAD INDEX ON A%02dB%02d AT (%+04d,%+04d) WHERE SHIFT [%+04d][%+04d] CAUSES (%+04d,%+04d)\n", label_a, label_b, i,j, shift_bx, shift_by, i+shift_bx, j+shift_by);
								printf(".	=>TRACE : SHIFT VALUE = [%+04d][%+04d]\n", shift_bx, shift_by);
								printf(".	.	=>TRACE : AT A%02dB%02d\n",label_a ,label_b);
								printf(".	.	.	=>TRACE : CENTRE VALUE = A%02d[%+04d][%+04d] B%02d[%+04d][%+04d]\n\n",label_b, centre_bx, centre_by, label_a, centre_ax, centre_ay);
							}

						}
					}
					else
					{
						//OUT OF BOUND ERROR
					}
				}
			}
			/*
			else if((shift_bx < shift_ax) && (shift_by < shift_ay))
			{
				if(image_label_b[i][j] == (label_a+L_BASE))
				{
					if(image_label_a[i+shift_ax][j+shift_ay] == (label_b+L_BASE))
					{
						total++;
						if((i+shift_ax < 0)||(j+shift_ay < 0))
						{
							printf("SHIFT ERROR!:BAD INDEX ON A%02dB%02d AT (%+04d,%+04d) WHERE SHIFT [%+04d][%+04d] CAUSES (%+04d,%+04d)\n", label_a, label_b, i,j, shift_ax, shift_ay, i+shift_ax, j+shift_ay);
							printf(".	=>TRACE : SHIFT VALUE = [%+04d][%+04d]\n", shift_ax, shift_ay);
							printf(".	.	=>TRACE : AT A%02dB%02d\n",label_a ,label_b);
							printf(".	.	.	=>TRACE : CENTRE VALUE = A%02d[%+04d][%+04d] B%02d[%+04d][%+04d]\n\n",label_b, centre_bx, centre_by, label_a, centre_ax, centre_ay);
						}
					}
				}
			}
			else if((shift_bx > shift_ax) && (shift_by < shift_ay))
			{
				if(image_label_b[i+shift_ax][j] == (label_a+L_BASE))
				{
					if(image_label_a[i][j+shift_ay] == (label_b+L_BASE))
					{
						total++;
						if((i+shift_ax < 0)||(j+shift_ay < 0))
						{
							printf("SHIFT ERROR!:BAD INDEX ON A%02dB%02d AT (%+04d,%+04d) WHERE SHIFT [%+04d][%+04d] CAUSES (%+04d,%+04d)\n", label_a, label_b, i,j, shift_ax, shift_ay, i+shift_ax, j+shift_ay);
							printf(".	=>TRACE : SHIFT VALUE = [%+04d][%+04d]\n", shift_ax, shift_ay);
							printf(".	.	=>TRACE : AT A%02dB%02d\n",label_a ,label_b);
							printf(".	.	.	=>TRACE : CENTRE VALUE = A%02d[%+04d][%+04d] B%02d[%+04d][%+04d]\n\n",label_b, centre_bx, centre_by, label_a, centre_ax, centre_ay);
						}
					}
				}
			}
			else if((shift_bx < shift_ax) && (shift_by > shift_ay))
			{
				if(image_label_b[i+shift_bx][j] == (label_a+L_BASE))
				{
					if(image_label_a[i][j+shift_ay] == (label_b+L_BASE))
					{
						total++;
						if((i+shift_bx < 0)||(j+shift_ay < 0))
						{
							printf("SHIFT ERROR!:BAD INDEX ON A%02dB%02d AT (%+04d,%+04d) WHERE SHIFT [%+04d][%+04d] CAUSES (%+04d,%+04d)\n", label_a, label_b, i,j, shift_ax, shift_ay, i+shift_ax, j+shift_ay);
							printf(".	=>TRACE : SHIFT VALUE = [%+04d][%+04d]\n", shift_ax, shift_ay);
							printf(".	.	=>TRACE : AT A%02dB%02d\n",label_a ,label_b);
							printf(".	.	.	=>TRACE : CENTRE VALUE = A%02d[%+04d][%+04d] B%02d[%+04d][%+04d]\n\n",label_b, centre_bx, centre_by, label_a, centre_ax, centre_ay);
						}
					}
				}
			}
			*/
			else
			{
				//printf("SHIFT ERROR!:BAD INDEX ON A%02dB%02d AT (%+04d,%+04d) WHERE SHIFT [%+04d][%+04d] CAUSES (%+04d,%+04d)\n", label_a, label_b, i,j, shift_ax, shift_ay, i+shift_ax, j+shift_ay);
				//printf(".	=>TRACE : SHIFT VALUE = [%+04d][%+04d]\n", shift_ax, shift_ay);
				//printf(".	.	=>TRACE : AT A%02dB%02d\n",label_a ,label_b);
				//printf(".	.	.	=>TRACE : CENTRE VALUE = A%02d[%+04d][%+04d] B%02d[%+04d][%+04d]\n\n",label_b, centre_bx, centre_by, label_a, centre_ax, centre_ay);

			}
			//MATRIX SHIFT
			//printf("SHIFT DATA: A%d => ASSIGNED INDEX => [%d][%d]\n", label_a, i, j);
			//printf("SHIFT DATA: B%d => ASSIGNED INDEX => [%d][%d]\n", label_b, i+shift_by, j+shift_bx);
			/*
			if((shift_bx > shift_ax)||(shift_by > shift_ay))
			{
				if(//image_label_a[i+shift_ax][j+shift_ax] == (label_a+L_BASE)
				//&&
				image_label_b[i][j] == (label_b+L_BASE)
				)
				{
					if(image_label_a[i+shift_ax][j+shift_ax] == (label_a+L_BASE))
					{
						total++;
						if((i+shift_ay < 0)||(j+shift_ax < 0))
							printf("SHIFT ERROR 01: A%d, B%d => BAD INDEX => %d,%d[%d][%d]\n", label_a, label_b, i,j, (shift_by), (shift_bx));
					}
				}

			}
			if((shift_bx < shift_ax)||(shift_by < shift_ay))//refinements here
			{
				//(i+shift_ay < 0)||(j+shift_ax < 0)
				if(//image_label_b[i+shift_by][j+shift_bx] == (label_b+L_BASE)
				//&&
				image_label_a[i][j] == (label_a+L_BASE)
				)
				{
					if(image_label_b[i+shift_by][j+shift_bx] == (label_b+L_BASE))
					{
						total++;
						if((i+shift_by < 0)||(j+shift_bx < 0))
							printf("SHIFT ERROR 02: A%d, B%d => BAD INDEX => %d,%d[%d][%d]\n", label_a, label_b, i,j, (shift_by), (shift_bx));
					}
				}

			}
			*/
			/*
			if(
					//A REMAINS STILL, NO CHANGE
					image_label_a[i][j] == (label_a+L_BASE) //&&
					//SHIFT B TO A, SHIFT PIXEL LOCATION

					//TODO ENSURE VALUE DOES NOT POINT TO UNDEFINED SEGMENTS!
					//TODO CHECK FOR NEGATIVE VALUE
					//TODO IF ELSE OR MIN MAX TO COMPENSATE NEGATIVE VALUES
			)
			{
			//INCREMENTAL
				if(image_label_b[i+shift_by][j+shift_bx] == (label_b+L_BASE))
				{
					total++;
					if((i+shift_by < 0)||(j+shift_bx < 0))
					{
						printf("SHIFT ERROR: A%d, B%d => ERROR INDEX => %d,%d[%d][%d]\n", label_a, label_b, i,j, (shift_by), (shift_bx));
					}
				}

			}
			*/
			//TODO ADD OTHER FUNCTIONS....
		}
	}
	//write_bmp_mono(image_label_and, "kirino_intersect.bmp");
	//write_bmp_mono(image_buf_a, "kirino_intersecta.bmp");
	//write_bmp_mono(image_buf_b, "kirino_intersectb.bmp");
	/*
	printf("[]=>TOTAL AB => %d; ",total);
	printf("SIZE A => %d; ",size_a);
	printf("SIZE B => %d; ",size_b);
	printf("LOGIC AB => %d/(%d + %d)-%d;",total, size_a, size_b, total);
	printf("MATH AB => %d/%d; \n",total, ((size_a+size_b)-total));
	*/
	//LOOP MATRIX CONTENT END
	double rtn = (double)total/(((double)size_a+(double)size_b)-(double)total);

	return rtn;

}



double features_moment(unsigned char image_in[Y_SIZE][X_SIZE], int label)
{
	int i, j;
	double moment = 0;
	for (i = 0; i < Y_SIZE; i++)
		{
			for (j = 0; j < X_SIZE; j++)
			{
				if(label != 0)
				{
					if(image_in[i][j] == label)
						moment = i*j*image_in[i][j];//CHECK LABLE
				}
				else
				{
					moment = i*j*image_in[i][j];
					// x*y*I(x,y)
				}

			}
		}
return moment;
}
void features_structure(unsigned char label_in[Y_SIZE][X_SIZE], double distance[], int cnt, int c_x[], int c_y[])
{
	int i, j;
	int c;
	//double d;




	for (i = 0; i < Y_SIZE; i++)
		{
			for (j = 0; j < X_SIZE; j++)
			{
				for(c = 0; c < cnt-1; c++)
				{
					distance[c] = sqrt(pow((c_x[c]-c_y[c]), 2) + pow((c_x[c+1]-c_y[c+1]), 2));

					//TODO DISTANCE FORMULA HERE....
					//TODO AND STRUCTURE STUFF....BLA BLA BLA
					//TODO CHECK K-MEANS.....
					//TODO GOOGLE STUFF ABOUT IMAGE IDENTIFICATION BLA BLA BLA

				}

			}
		}
}

double general_segmentation(unsigned char image_a[3][Y_SIZE][X_SIZE],unsigned char image_b[3][Y_SIZE][X_SIZE], int segment_size, int mode, int border, double diffval[4][3])
{
	int i, j;
	int c;
	int array_size = segment_size*segment_size;

	div_t segment_y = div(Y_SIZE, segment_size);
	//segment_y;
	div_t segment_x = div(X_SIZE, segment_size);
	/*
    t_array_index ***position_neighbor;// =  malloc(sizeof ***position_neighbor * array_size * 3 * array_size);
    position_neighbor =  malloc(sizeof ***position_neighbor * array_size * 3 * array_size);
    for (c = 0; c < 3; c++)
    {
    	position_neighbor[c] =  malloc(sizeof(t_array_index) * array_size * array_size);
        for (i = 0; i < array_size; i++)
        {
        	position_neighbor[c][i] =  malloc(sizeof(t_array_index) * array_size);
        }
    }

    t_array_index ***neighbor_coefficient_deviation;// = neighbor_coefficient_deviation = malloc(sizeof ***neighbor_coefficient_deviation * array_size * 3 * array_size);
    neighbor_coefficient_deviation = malloc(sizeof ***neighbor_coefficient_deviation * array_size * 3 * array_size);
    for (c = 0; c < 3; c++)
    {
    	neighbor_coefficient_deviation[c] = malloc(sizeof(t_array_index) * array_size * array_size);
        for (i = 0; i < array_size; i++)
        {
        	neighbor_coefficient_deviation[c][i] = malloc(sizeof(t_array_index) * array_size);
        }
    }
	 */



	//double d;

	//segment_x;
	static int matrix_a[3];
	static int matrix_b[3];
	static double diff;
	static double mn_diff[3];
	static double mx_diff[3];
	static double min_diff;
	static double max_diff;

	unsigned char image_out[3][Y_SIZE][X_SIZE];
	//int segment_area = ((segment_x.quot)*(segment_y.quot));
	double matrix_diff[3][array_size][2];
	//double matrix_diff[3][segment_size*segment_size][2];
	double seglength_x = segment_x.quot;
	double seglength_y = segment_y.quot;

	double search_matrix[3][array_size][array_size];
	double link_matrix[2][3][array_size][array_size];

	//double search_identifier[3][1][2];





	double ratio_matrix[3][array_size][array_size];
	double ratio_identifier[3][array_size];
	double ratio_unity[3][array_size];
	static double ratio_deviation[3];
	static double ratio_average[3];
	static double channel_diff[3];

	int structure[2][array_size];
	//int array_size = segment_size*segment_size;
    int neighbor_range = 8;

	if((border*2 >= seglength_x) || (border*2 >= seglength_y))
		border = roundf(((float)(fmin(seglength_x,seglength_y)-1))/2);
	double neighbor_coefficient_matrix[3][array_size][array_size][neighbor_range+1];



	//advance position_neighbor part




	t_array_index position_neighbor[3][array_size][array_size];
	t_array_index_ext neighbor_coefficient_deviation[3][array_size][array_size];
    //double neighbor_coefficient_deviation[array_size][array_size];
    double segment_coefficient_deviation[3][array_size];
    double segment_coefficient_level[c];

    double structure_deviation[3][array_size/2][array_size];
    double structure_coefficient[3][array_size/2];
    //double neighbor_coefficient_average[segment_size*segment_size][segment_size*segment_size][neighbor_range];


/*	int neighbor_dimension[] = {-1, 1,

			(-(segment_size)), (segment_size),

			(-(segment_size-1)), (segment_size+1),

			(-(segment_size+1)), (segment_size-1),

			(-(segment_size*2)), (segment_size*2),

			-2, 2};
*/
	int neighbor_dimension[] = {-1, 1,

			(-(segment_size)), (segment_size),

			(-(segment_size-1)), (segment_size+1),

			(-(segment_size+1)), (segment_size-1)};
	//could be (n=1)8, (n=2)24-> 8+(8*2)....or even a circle...
	//this one is +1 expansion in x and y dimension only,thus, 4



    //



	//int segemnt = Y_SIZE/segment_size;
	int seg_y, seg_x;
	int counter = 0;
	int pixel_c = 0;
	int pixel_s = 0;

	for(seg_y = 0; seg_y < segment_size; seg_y++)
	{
		for(seg_x = 0; seg_x < segment_size; seg_x++)
		{
			for (i = (seglength_y*seg_y)+border; i < ((seglength_y)*(seg_y+1))-border; i++)
			{
				for (j = (seglength_x*seg_x)+border; j < ((seglength_x)*(seg_x+1))-border; j++)
				{
					for(c = 0;c < 3;c++)
					{
						matrix_a[c] += image_a[c][i][j];
						matrix_b[c] += image_b[c][i][j];
					}
					pixel_c++;
					pixel_s++;
				}
			}
			for(c = 0;c < 3;c++)
			{
				matrix_diff[c][counter][0] = (double)matrix_a[c]/((seglength_x-(border*2))*(seglength_y-(border*2)));
				matrix_diff[c][counter][1] = (double)matrix_b[c]/((seglength_x-(border*2))*(seglength_y-(border*2)));
				matrix_a[c] = 0;
				matrix_b[c] = 0;
			}
			counter++;
			if(pixel_s == ((seglength_x-(border*2))*(seglength_y-(border*2))))
			{
				pixel_s = 0;
			}
			else
			{
				//printf("BAD INDEX AT %d WHERE [%d] != [%f]\n",counter, pixel_s,((seglength_x-(border*2))*(seglength_y-(border*2))));
				pixel_s = 0;
			}
		}
	}

	double diff_avg = (double)0;
	if(mode == 1)
	{
		for(c = 0;c < 3;c++)
		{
			channel_diff[c] = euclidean_dist(matrix_diff[c], array_size);
		}
		diff_avg = stat_stddev(channel_diff, 3);
	}

	for(c = 0;c < 3;c++)
	{
		diff += euclidean_dist(matrix_diff[c], array_size);
	}
	//section above is border independent
	//section below is area and pixel dependent(aka border dependent)
	int du, di;
	for(du = 0; du < segment_size*segment_size; du++)
	{
		//printf("[%d]=>{", du);
		for(di = 0; di < segment_size*segment_size; di++)
		{
			for(c = 0;c < 3;c++)
			{
				//search_identifier[c][0][0] = matrix_diff[c][du][0];
				//search_identifier[c][0][1] = matrix_diff[c][di][1];
				//search_matrix[c][du][di] = euclidean_dist(search_identifier[c], 1);
				//self relevance comparison, ex: |A-B|, |A-C| .....
				search_matrix[c][du][di] = fabs(matrix_diff[c][du][0] - matrix_diff[c][di][1]);
				link_matrix[0][c][du][di] = fabs(matrix_diff[c][du][0] - matrix_diff[c][di][0]);
				link_matrix[1][c][du][di] = fabs(matrix_diff[c][du][1] - matrix_diff[c][di][1]);

				//TODO we need all channels here, otherwise, the only channel left would be blue since 0r1g2b
				position_neighbor[c][du][di].value = fabs(matrix_diff[c][du][0] - matrix_diff[c][di][1]);
				position_neighbor[c][du][di].index = di;

				//printf("C%d,A%d:B%d[%f];",c,du,di,search_matrix[c][du][di]);
			}
			//printf("\n");
		}
		//printf("}\n");
		for(c = 0;c < 3;c++)
			qsort(search_matrix[c][du], sizeof(search_matrix[c][du])/sizeof(search_matrix[c][du][0]), sizeof(search_matrix[c][du][0]), cmp);

		//sort array index(struct index)
		for(c = 0;c < 3;c++)
			qsort(position_neighbor[c][du], array_size, sizeof(position_neighbor[c][du][0]), cmp_struct);

	}
	/*
	for(du = 0; du < segment_size*segment_size; du++)
	{
		for(di = (segment_size*segment_size)-1; di >= 0; di--)
		{
			//neighbor_coefficient_deviation[du][di] = stat_stddev(neighbor_coefficient_matrix[du][di], neighbor_range);
			//printf("LIST ARRAY -> %d => %f\n", position_neighbor[du][di].index, position_neighbor[du][di].value);
		}
	}
	*/

	counter = 0;
	//advance search!
	/*
	for(du = 0; du < segment_size*segment_size; du++)
	{
		//matrix_diff[c][du][0]

	}
	*/

	c = 1;
	for(du = 0; du < segment_size*segment_size; du++)
	{
		//neighbor_matrix[du] = ;

		for(di = (segment_size*segment_size)-1; di >= 0; di--)
		{

			for(c = 0; c < 3; c++)
			{
				//search_matrix[c][du];
				//segment neighbor loop
				for(i = 0; i < neighbor_range; i++)
				{
					if((du+neighbor_dimension[i] > 0) && (position_neighbor[c][du][di].index + neighbor_dimension[i] > 0))
					{
						neighbor_coefficient_matrix[c][du][di][i] = fabs( matrix_diff[c][du+neighbor_dimension[i]][0] - matrix_diff[c][position_neighbor[c][du][di].index + neighbor_dimension[i]][1] );
						/*
						neighbor_coefficient[i] = fmin(matrix_diff[c][position_neighbor[du][di].index+neighbor_dimension[i]][0], matrix_diff[c][position_neighbor[du][di].index+neighbor_dimension[i]][1])/
								fmin(matrix_diff[c][position_neighbor[du][di].index+neighbor_dimension[i]][0], matrix_diff[c][position_neighbor[du][di].index+neighbor_dimension[i]][1]);
						 */
						//printf("[%d] -> abs(%f-%f)=%f; POS00:%d POS01:%d \n", i, matrix_diff[c][du+neighbor_dimension[i]][0], matrix_diff[c][position_neighbor[c][du][di].index + neighbor_dimension[i]][1],
						//	neighbor_coefficient_matrix[c][du][di][i], du+neighbor_dimension[i], position_neighbor[c][du][di].index + neighbor_dimension[i] );
						counter++;
					}
				}
				counter++;
				neighbor_coefficient_matrix[c][du][di][neighbor_range] = position_neighbor[c][du][di].value;

				//printf("[%d] -> abs(%f-%f)=%f; POS00:%d POS01:%d \n", 12, matrix_diff[c][du][0], matrix_diff[c][position_neighbor[c][du][di].index][1],
				//	neighbor_coefficient_matrix[c][du][di][neighbor_range], du, position_neighbor[c][du][di].index);

				neighbor_coefficient_deviation[c][du][di].value = stat_stddev(neighbor_coefficient_matrix[c][du][di], counter);

				//neighbor_coefficient_deviation[du][di].value = position_neighbor[du][di].index;
				neighbor_coefficient_deviation[c][du][di].index = du;
				neighbor_coefficient_deviation[c][du][di].index_assoc = position_neighbor[c][du][di].index;
				//if(neighbor_coefficient_deviation[du][(segment_size*segment_size)-1] == (double)0)
					//printf("ADVANCE SEARCH:[%d(%d)][%d] -> %f\n", du, position_neighbor[du][di].index+neighbor_dimension[i], di, neighbor_coefficient_deviation[du][di]);

				//qsort(position_neighbor[du], array_size, sizeof(position_neighbor[du][0]), cmp_struct);

				//segment_coefficient_deviation[du] = neighbor_coefficient_deviation[du][(segment_size*segment_size)-1].value;
				/*
				if(neighbor_coefficient_deviation[c][du][di].value == 0)
					printf("ADVANCE SEARCH:[%d(%d)][%d] -> %f << peak\n", du, position_neighbor[c][du][di].index, di, neighbor_coefficient_deviation[c][du][di].value);
				if(neighbor_coefficient_deviation[c][du][di].value != 0)
					printf("ADVANCE SEARCH:[%d(%d)][%d] -> %f\n", du, position_neighbor[c][du][di].index, di, neighbor_coefficient_deviation[c][du][di].value);
				 */
				if(neighbor_coefficient_deviation[c][du][di].value == 0)
				{
					//structure[0] = du;
					//structure[1] = position_neighbor[c][du][di].index;

				}

				//SEGMENT W.A.S
				counter = 0;
			}

		}
		//counter++;
		//TODO SORT NEIGHBOR ARRAY HERE WITH QSORT: THE ORIGINAL DATA FROM ONE TO ONE MIGHT LOST ACCURACY...
		for(c = 0; c < 3; c++)
			qsort(neighbor_coefficient_deviation[c][du], array_size, sizeof(neighbor_coefficient_deviation[c][du][0]), cmp_struct);
		for(i = 0; i < neighbor_range; i++)
		{

		}
		for(c = 0; c < 3; c++)
			//printf("COEFFICIENT PEAK:%f\n",neighbor_coefficient_deviation[c][du][(segment_size*segment_size)-1].value);
		for(c = 0; c < 3; c++)
			segment_coefficient_deviation[c][du] = neighbor_coefficient_deviation[c][du][array_size-1].value;
		for(c = 0; c < 3; c++)
			segment_coefficient_level[c] += neighbor_coefficient_deviation[c][du][array_size-1].value;
		for(c = 0; c < 3; c++)
		{
			for(i = 0; i < 8; i++)
			{
				structure_deviation[c][i][du] = fabs(neighbor_coefficient_deviation[c][du][array_size-i].index-neighbor_coefficient_deviation[c][du][array_size-i].index_assoc);
				//structure_level[c][i][du] = fabs(neighbor_coefficient_deviation[c][du][array_size-i].index-neighbor_coefficient_deviation[c][du][array_size-i].index_assoc);
				//printf("AB LOCK: %f\n", structure_deviation[c][i][du]);

			}
		}

	}

	for(c = 0; c < 3; c++)
	{
		for(i = 0; i < 8; i++)
		{
			structure_coefficient[c][i] = stat_stddev(structure_deviation[c][i], array_size);
		}
		qsort(structure_coefficient[c], sizeof(structure_coefficient[c])/sizeof(structure_coefficient[c][0]), sizeof(structure_coefficient[c][0]), cmp);
		//printf("STRUCT LOCK: %f\n", structure_coefficient[c][0]);
	}




	//printf("SEGMENT COEFFICIENT:%f\n",stat_stddev(segment_coefficient_deviation, segment_size*segment_size));










	//advance search!
	counter = 0;
	for(du = 0; du < segment_size*segment_size; du++)
	{
		for(di = 0; di < segment_size*segment_size; di++)
		{
			for(c = 0;c < 3;c++)
			{
				//self & other ratio......
				//link_matrix[0][c][du][di];
				//link_matrix[1][c][du][di];
				if((link_matrix[0][c][du][di] != (double)0) || (link_matrix[0][c][du][di] != (double)0))
				{
					//printf("C%d,Ri:%f; =>[%f][%f] \n",c, ratio_matrix[c][du][di],link_matrix[0][c][du][di],link_matrix[1][c][du][di]);
					ratio_matrix[c][du][di] = (fmin(link_matrix[0][c][du][di], link_matrix[1][c][du][di])/fmax(link_matrix[0][c][du][di], link_matrix[1][c][du][di]));
					//if(du == (double)0)
						ratio_identifier[c][du] = (double)1;

					ratio_identifier[c][du]  *= (fmin(link_matrix[0][c][du][di], link_matrix[1][c][du][di])/fmax(link_matrix[0][c][du][di], link_matrix[1][c][du][di]));
					counter++;
					//printf("C%d;du[%d] =>r[%f] \n",c, du,(fmin(link_matrix[0][c][du][di], link_matrix[1][c][du][di])/fmax(link_matrix[0][c][du][di], link_matrix[1][c][du][di])));
				}
			}
		}
		for(c = 0;c < 3;c++)
		{
			//ratio_unity[c][du]  = ratio_identifier[c][du]/((segment_size*segment_size)-1);
			//if(ratio_identifier[c][du] == (double)1)
				//ratio_identifier[c][du] = 0;
			ratio_unity[c][du]  = ratio_identifier[c][du];
			//printf("!C%d;du[%d] =>r[%f] \n",c, du,ratio_identifier[c][du]);
			//printf("C%d,Ru:%f;\n",c, ratio_identifier[c][du]);
			//printf("C%d,Ra:%f;\n",c, ratio_unity[c][du]);
			//ratio_deviation[c][du] /= (segment_size*segment_size);
			//ratio_deviation[c][du]

		}
	}
	for(c = 0;c < 3;c++)
	{
		ratio_deviation[c] = stat_stddev(ratio_unity[c], (segment_size*segment_size));
		//printf("C%d,Rd:%f,",c, ratio_deviation[c]);
		ratio_average[c] = stat_avg(ratio_unity[c], (segment_size*segment_size));
		//printf("Ra:%f;\n", ratio_average[c]);
		//printf("Ro:%f;\n", ratio_average[c]);
	}
	//printf("Count:%d;\n",counter);





	for(i = 0; i < segment_size*segment_size; i++)
	{
		//printf("[%d] => ", i);
		for(c = 0;c < 3;c++)
		{
			//printf("[%d][%f];\n", c, search_matrix[c][i][(segment_size*segment_size)-1]);
			mn_diff[c] += search_matrix[c][i][(segment_size*segment_size)-1];
			mx_diff[c] += search_matrix[c][i][0];

		}
		//printf("\n");

	}

	for(c = 0;c < 3;c++)
	{
		//printf("[%f];\n", mn_diff[c]);
		//printf("[%f];\n", mx_diff[c]);
		min_diff += (mn_diff[c]/3);
		max_diff += (mx_diff[c]/3);
	}


	/*
	printf("\nEUCLIDEAN DISTANCE:\n");
	printf("SEGMENTATION ASSIGNED   :%d\n",segment_size*segment_size);
	printf("SEGMENTATION TOTAL      :%d\n",counter);
	printf("SEGMENT PIXEL ASSIGNED  :%f\n",(segment_size*seglength_x)*(segment_size*seglength_y));
	printf("SEGMENT PIXEL TOTAL     :%d\n",pixel_c);

	printf("SEGMENTATION LENGTH X   :%f\n",seglength_x);
	printf("SEGMENTATION LENGTH Y   :%f\n",seglength_y);
	printf("SEGMENTATION AREA       :%f\n",(seglength_x*seglength_y));
	printf("R                       :%f\n",euclidean_dist(matrix_diff[0], segment_size*segment_size));
	printf("G                       :%f\n",euclidean_dist(matrix_diff[1], segment_size*segment_size));
	printf("B                       :%f\n",euclidean_dist(matrix_diff[2], segment_size*segment_size));
	printf("CHANNEL DEVIATION       :%f\n",diff_avg);

	printf("MIN POSSIBLE DISTANCE   :%f\n", min_diff/(segment_size*segment_size));
	printf("MAX POSSIBLE DISTANCE   :%f\n", max_diff/(segment_size*segment_size));
	printf("BORDER                  :%d\n", border);
	*/


	counter = 0;
	double matrix_identifier[1][2];

	for(seg_y = 0; seg_y < segment_size; seg_y++)
	{
		for(seg_x = 0; seg_x < segment_size; seg_x++)
		{
			for (i = seglength_y*seg_y; i < ((seglength_y)*(seg_y+1)); i++)
			{
				for (j = seglength_x*seg_x; j < ((seglength_x)*(seg_x+1)); j++)
				{
					for(c = 0;c < 3;c++)
					{
						//image_out[c][i][j] = segment_identifier;

						matrix_identifier[0][0] = matrix_diff[c][counter][0];
						matrix_identifier[0][1] = matrix_diff[c][counter][1];
						//image_out[c][i][j] = abs(image_a[c][i][j] - euclidean_dist(matrix_identifier, 1));
						image_out[c][i][j] = matrix_diff[c][counter][0];
					}
				}
			}
			counter++;
		}
	}


	//qsort(d[a], cnt_b, sizeof(double), sort);
	//write_bmp_color(image_out, "SEGMENTATION_VISUAL.bmp");
	//diff = euclidean_dist(matrix_diff, segment_size*segment_size);
	diffval[0][0] = min_diff/(segment_size*segment_size);
	diffval[0][1] = max_diff/(segment_size*segment_size);
	diffval[0][2] = diff/3;

	for(c = 0;c < 3;c++)
		diffval[1][c] = ratio_average[c];
	for(c = 0;c < 3;c++)
		diffval[2][c] = ratio_deviation[c];
	for(c = 0;c < 3;c++)
		diffval[3][c] = segment_coefficient_level[c]*(structure_coefficient[c][7]);

	//diffval[3][c] = stat_stddev(segment_coefficient_deviation[c], segment_size*segment_size);


	//free(position_neighbor);
	//free(neighbor_coefficient_deviation);

	//diffval[1][1] = ratio_deviation;
	return 0;
}
/*
double general_distance(unsigned char image_a[3][Y_SIZE][X_SIZE],unsigned char image_b[3][Y_SIZE][X_SIZE])
{
	int i, j;
	int c;
	//int matrix_a[3],matrix_b[3];
	double diff;
	double matrix_diff[3][Y_SIZE*X_SIZE][2];
	//int seg_y, seg_x;
	int counter;
	for (i = 0; i < Y_SIZE; i++)
	{
		for (j = 0; j < Y_SIZE; j++)
		{
			for(c = 0;c < 3;c++)
			{
				matrix_diff[c][counter][0] = image_a[c][i][j];
				matrix_diff[c][counter][1] = image_b[c][i][j];
			}
			counter++;
		}
	}
	for(c = 0;c < 3;c++)
	{
		diff += euclidean_dist(matrix_diff[c], counter);
	}
	printf("\nEUCLIDEAN DISTANCE:\n");
	printf("SEGMENTATION TOTAL:%d\n",counter);
	printf("R:%f\n",euclidean_dist(matrix_diff[0], counter));
	printf("G:%f\n",euclidean_dist(matrix_diff[1], counter));
	printf("B:%f\n",euclidean_dist(matrix_diff[2], counter));
	//diff = euclidean_dist(matrix_diff, segment_size*segment_size);
	return diff/3;
}
*/
void print_image(unsigned char image[3][Y_SIZE][X_SIZE],int channel)
{
	int i, j;
	for (i = 0; i < Y_SIZE; i++)
	{
		printf("\n");
		for (j = 0; j < Y_SIZE; j++)
		{
			printf("[%d]",image[channel][i][j]);
		}
	}
}
/*
double calc_equalize(unsigned char image_label[Y_SIZE][X_SIZE],
		unsigned char image_label_equalized[Y_SIZE][X_SIZE],
	int label, int cnt_o, int cnt_eq)
{
	int		i, j;
	int cnt;
	double	 total;
	int size[HIGH];
	//i+L_BASE
	total = 0;
	for(cnt = 0; cnt < cnt_o; cnt++)
	{
		for (i = 0; i < Y_SIZE; i++)
		{
			for (j = 0; j < X_SIZE; j++)
			{
				if (image_label[i][j] == cnt+L_BASE)
				{
					//DEFINE L_BASE as BASE SERIAL
					total++;
					//CALC SIZE
				}
			}
		}
		size[cnt] = total;
	}

	if (total == 0.0) return 0.0;

	return total;
}
*/
/*--- calc_length ----------*/
double calc_length(unsigned char image_label[Y_SIZE][X_SIZE], int label)
{
	int		i, j;
	double	trace();

	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			if (image_label[i][j] == label)	return trace(image_label, j-1, i);
	return 0;
}

/*--- trace ---------*/
double trace(unsigned char image_label[Y_SIZE][X_SIZE], int xs, int ys)
{
	int		x, y, no, vec;
	double	l;

	l = 0;	x = xs; y = ys; no = image_label[y][x+1]; vec = 5;
	for (;;) {
		if (x == xs && y == ys && l != 0) return l;
		image_label[y][x] = HIGH;
		switch (vec) {
			case 3:
				if (image_label[y][x+1] != no && image_label[y-1][x+1] == no)
					{x = x+1;   ; l++       ; vec = 0; continue;}
					/* no break */
			case 4:
				if (image_label[y-1][x+1] != no && image_label[y-1][x] == no)
					{x = x+1; y = y-1; l += ROOT2; vec = 1; continue;}
					/* no break */
			case 5:
				if (image_label[y-1][x] != no && image_label[y-1][x-1] == no)
					{  ; y = y-1; l++       ; vec = 2; continue;}
					/* no break */
			case 6:
				if (image_label[y-1][x-1] != no && image_label[y][x-1] == no)
					{x = x-1; y = y-1; l += ROOT2; vec = 3; continue;}
					/* no break */
			case 7:
				if (image_label[y][x-1] != no && image_label[y+1][x-1] == no)
					{x = x-1;   ; l++       ; vec = 4; continue;}
					/* no break */
			case 0:
				if (image_label[y+1][x-1] != no && image_label[y+1][x] == no)
					{x = x-1; y = y+1; l += ROOT2; vec = 5; continue;}
					/* no break */
			case 1:
				if (image_label[y+1][x] != no && image_label[y+1][x+1] == no)
					{  ; y = y+1; l++       ; vec = 6; continue;}
					/* no break */
			case 2:
				if (image_label[y+1][x+1] != no && image_label[y][x+1] == no)
					{x = x+1; y = y+1; l += ROOT2; vec = 7; continue;}
				vec = 3;
				/* no break */
		}
	}
	return 0;
}




/*--- extract_ratio --------*/
void extract_ratio(unsigned char image_label_in[Y_SIZE][X_SIZE],
	unsigned char image_label_out[Y_SIZE][X_SIZE],
	int cnt, double ratio[], double ratio_min, double ratio_max)
{
	int	i, j, x, y;
	int	lno[256];

	for (i = 0, j = 0; i < cnt; i++)
	    if (ratio[i] >= ratio_min && ratio[i] <= ratio_max)
			lno[j++] = L_BASE+i;
	for (y = 0 ; y < Y_SIZE; y++) {
		for (x = 0; x < X_SIZE; x++) {
			image_label_out[y][x] = 0;
			for (i = 0; i < j; i++)
			    if (image_label_in[y][x] == lno[i])
			    	image_label_out[y][x] = image_label_in[y][x];
	    }
	}
}


/*---- masking --- ----*/
void masking(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE],
	unsigned char image_mask[Y_SIZE][X_SIZE])
{
	int	i, j;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
		  	if (image_mask[i][j] == HIGH)	image_out[i][j] = image_in[i][j];
		  	else							image_out[i][j] = 0;
		}
	}
}
void label_masking(unsigned char image_lable[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE],
	int label[], int cc)
{
	int	i, j;
	int c;
	//printf("LABEL => \n");
	for(c = 1; c < cc+1; c++)
		//printf("MAKING LABEL => %d\n", label[c]);

	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			image_out[i][j] = 0;

	for (i = 0; i < Y_SIZE; i++)
	{
		for (j = 0; j < X_SIZE; j++)
		{
			for(c = 1; c < cc+1; c++)
			{
				//printf("MAKING LABEL => %d", label[c+L_BASE]);
				if (label[c] == image_lable[i][j])	image_out[i][j] = HIGH;
				//else							image_out[i][j] = 0;
			}

		}
	}
}





/*--- extract_size --------*/
void extract_size(unsigned char image_label_in[Y_SIZE][X_SIZE],
	unsigned char image_label_out[Y_SIZE][X_SIZE],
	int cnt, double size[], double size_min, double size_max)
{
	int	i, j, x, y;
	int	lno[256];

	for (i = 0, j = 0; i < cnt; i++)
		if (size[i] >= size_min && size[i] <= size_max)	lno[j++] = L_BASE+i;
	for (y = 0; y < Y_SIZE; y++) {
	    for (x = 0; x < X_SIZE; x++) {
			image_label_out[y][x] = 0;
	    	for (i=0 ; i<j ; i++)
		    	if (image_label_in[y][x] == lno[i])
		    		image_label_out[y][x] = image_label_in[y][x];
		}
	}
}

/*
 * FEATURE END
 */

/*
 * NOISE
 */


/*--- median -------------------------------*/
void median(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE])
{
	int i, j;
  	unsigned char c[9];

	for (i = 1; i < Y_SIZE-1; i++) {
		for (j = 1; j < X_SIZE-1; j++) {
			c[0] = image_in[i-1][j-1];
			c[1] = image_in[i-1][j];
			c[2] = image_in[i-1][j+1];
			c[3] = image_in[i][j-1];
			c[4] = image_in[i][j];
			c[5] = image_in[i][j+1];
			c[6] = image_in[i+1][j-1];
			c[7] = image_in[i+1][j];
			c[8] = image_in[i+1][j+1];
			image_out[i][j] = median_value(c);
		}
	}
}
/*--- median_value -------------------------------------------*/
int median_value(unsigned char c[9])
{
	int i, j, buf;

	for (j = 0; j < 8; j++) {
		for (i = 0; i < 8; i++) {
			if (c[i+1] < c[i]) {
				buf = c[i+1];
				c[i+1] = c[i];
				c[i] = buf;
			}
		}
	}
	return c[4];
}

void image_multiplication_color(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	unsigned char image_out[3][Y_SIZE][X_SIZE])
{
	int	i, j, k;
	double d;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < Y_SIZE; j++) {
			for (k = 0; k < X_SIZE; k++) {
				d = (double)image_in1[i][j][k] * image_in2[i][j][k] / 255.0;
				image_out[i][j][k] = (unsigned char)d;
			}
		}
	}
}


/*
 * NOISE END
 */
/*
 * CORE UTILITY
 */
/*
void core_trace(const char* d1, const char* d2)
{

}
void core_log(const int level , const char* data)
{

}
*/

/*
 * CORE UTILITY END
 */

