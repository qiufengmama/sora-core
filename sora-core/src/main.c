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
	unsigned char image_label[2][Y_SIZE][X_SIZE];
	unsigned char image_feature[2][Y_SIZE][X_SIZE];
	unsigned char image_sub[2][3][Y_SIZE][X_SIZE];
	unsigned char image_out[3][Y_SIZE][X_SIZE];
	unsigned char image_test[2][3][Y_SIZE][X_SIZE];
	//char *image_path = "miyako.bmp";
	//char *image_0_path = "miyako_FFT.bmp";
	//char *image_1_path = "miyako_FFT_flt.bmp";
	//image_path = "MIYAKO.bmp";
	//clock_t start = clock();
	//clock_t start_read = clock();
	//clock_t start_fft = clock();
/*
 * VAR DEFINE
 */
	//LISTING DEFAULT
	int threshold_rgb = 50;
	int threshold_ysh = 50;
	double threshold_ff = 0.5; // *100
	int lim = (Y_SIZE+X_SIZE)/2;
	int a = 0;
	int b = 128;
	int clrv = 1;
	//DEBUG ONLY

	//debug(str);


	char *path_fft_a = join(join(join(CORE_CODE_NAME, "_FFT_A_"), tochar(TIME)), ".bmp");
	char *path_fft_b = join(join(join(CORE_CODE_NAME, "_FFT_B_"), tochar(TIME)), ".bmp");
	char *path_fft_x = join(join(join(CORE_CODE_NAME, "_FFT_X_"), tochar(TIME)), ".bmp");

	char *path_ff_a = join(join(join(CORE_CODE_NAME, "_FF_A_"), tochar(TIME)), ".bmp");
	char *path_ff_b = join(join(join(CORE_CODE_NAME, "_FF_B_"), tochar(TIME)), ".bmp");
	char *path_ff_x = join(join(join(CORE_CODE_NAME, "_FF_X_"), tochar(TIME)), ".bmp");

	char *path_thr_a = join(join(join(CORE_CODE_NAME, "_THR_A_"), tochar(TIME)), ".bmp");
	char *path_thr_b = join(join(join(CORE_CODE_NAME, "_THR_B_"), tochar(TIME)), ".bmp");

	int debug = 1;
	char *debug_a = "ko.bmp";
	char *debug_b = "miyako.bmp";
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
	//clock_t start_arr = clock();

/*
 * DATA READ
 */
	//DECL
if(debug == 0)
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
	long fslen = (((lim*lim)*3)*3+((lim*lim)*3))*2+lim;
	char c[fslen];
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
}
else if(debug  == 1)
{
	//image_test[ii][ci][xi][yi]
	read_bmp_color(image_test[0], debug_a); // <-
	read_bmp_color(image_test[1], debug_b); // <-
}
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
char text_buf_a[65536];
char text_buf_b[65536];
char buf_f[65536];
unsigned char image_label_out_a[Y_SIZE][X_SIZE];
unsigned char image_label_out_b[Y_SIZE][X_SIZE];
double size_a[HIGH], ratio_a[HIGH], length_a[HIGH];
double size_b[HIGH], ratio_b[HIGH], length_b[HIGH];
double size_f[HIGH], ratio_f[HIGH], length_f[HIGH];
unsigned char image_ffl[2][Y_SIZE][X_SIZE];
cnt_a = labeling(image_thresh[0], image_label[0]); //<<
cnt_b = labeling(image_thresh[1], image_label[1]); //<<
assign(image_ffl[0], image_label[0]);// <-
assign(image_ffl[1], image_label[1]);// <-
//::

//labeling(image_thresh[0], image_label[0], *cnt, *buf);
features(image_label[0], image_feature[0], cnt_a, size_a, length_a, ratio_a, text_buf_a); // ->.
features(image_label[1], image_feature[1], cnt_b, size_b, length_b, ratio_b, text_buf_b); // ->.

//printf("FEATURE A:%d\n", ;

//printf("CNT A:%d\n", cnt_a);
//printf("CNT B:%d\n", cnt_b);
//printf("\nA=>\n%s\n<=A\n", text_buf_a);
//printf("\nB=>\n%s\n<=B\n", text_buf_b);


double rslt;
int label[HIGH];

features_compare(image_ffl[0],
	image_ffl[1],
	image_label_out_a,
	image_label_out_b,
	cnt_a, cnt_b, size_f,
	length_f, ratio_f, buf_f, &rslt, label, &cc, threshold_ff); // ->.

label_masking(image_label[0], image_marker, label, cc); // ->

//printf("\nCC=>\n%s\n<=CC\n", buf_f);
//double *crt_rtn = &rslt;
//printf("=>LAB: %d \n", &label);
printf("%f;\n", rslt);


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
int cl,i,j,fl;
int true_val=0,false_val=0,total_val=0;

if(debug  == 1)
{
write_bmp_color(image_sub[0], path_fft_a);
write_bmp_color(image_sub[1], path_fft_b);
}



//if((WRITE_FFT0 != -1) && (WRITE_FFT1 != -1))
//{//check FFT
	//printf("FFT WRITE OK!\n");
	//printf("WRITE Time elapsed: %fs\n", ((double)clock() - start_write) / CLOCKS_PER_SEC);
	//printf("TOTAL Time elapsed: %fs\n", ((double)clock() - start) / CLOCKS_PER_SEC);}
	rgb_to_ysh(image_sub[0],image_ysh[0]);//ysh convert image 0
	rgb_to_ysh(image_sub[1],image_ysh[1]);//ysh convert image 1
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
				for(cl = 0;cl <=2; cl++)
				{
					//RGB layer
					//exclude 0(!filtered)
					//RGB DISTANCE
					image_dist[0][cl][i][j] = abs(image_sub[0][cl][i][j]-image_sub[1][cl][i][j]);
					//YSH DISTANCE
					image_dist[1][cl][i][j] = abs(image_ysh[0][cl][i][j]-image_ysh[1][cl][i][j]);
					//DISTANCE WRITE
				}
			//if((image_sub[0][1][i][j]  != 0)){
				if(
						/*(image_sub[0][1][i][j]  != 0) && */
						(((image_dist[0][0][i][j] +
						image_dist[0][1][i][j] +
						image_dist[0][2][i][j]) <= threshold_rgb) ||
						((image_dist[1][0][i][j] <= threshold_ysh) &&
						(image_dist[1][1][i][j] <= threshold_ysh) &&
						(image_dist[1][2][i][j] <= threshold_ysh)))
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
unsigned char image_conv[clrv+1][3][Y_SIZE][X_SIZE];
assign(image_conv[0][0], image_rslt[0]);
assign(image_conv[0][1], image_rslt[1]);
assign(image_conv[0][2], image_rslt[2]);
for(clr = 0;clr <= clrv;clr++)
{
	convolute(image_conv[clr][0], image_conv[clr+1][0], b);
	convolute(image_conv[clr][1], image_conv[clr+1][1], b);
	convolute(image_conv[clr][2], image_conv[clr+1][2], b);
}
assign(image_out[0], image_conv[clrv+1][0]);
assign(image_out[1], image_conv[clrv+1][1]);
assign(image_out[2], image_conv[clrv+1][2]);
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

double c_per = (double)c_true/c_total;
printf("%f", c_per);
/*
printf("PERCENTAGE SIMILAR: %f%%,\n", c_per*100);
printf("PIXEL TOTAL: %d,\n", c_total);
printf("PIXEL SIMILAR: %d,\n", c_true);
printf("PIXEL DIFFER: %d,\n", c_false);
*/

if(debug  == 1)
{
write_bmp_color(image_out, path_fft_x);
}


/*
//new blank end
	int arr_i;
	int abi;
	int ax; //y
	int bx; //x

	for(arr_i = 0;arr_i <= 2; arr_i++){
		//printf("ARRAY %d  (", arr_i);
	for(ax = 0; ax <= lim-1;ax++){
		for(bx = 0;bx <= lim-1;bx++){
			//if(image_test[abi][arr_i][ax][bx] != 0)
				printf("%d,", image_main[arr_i][ax][bx]);
			//puts(image_test[abi][arr_i][ax][bx]);
		}
	}
	//printf(");");
	printf("|\n");
	}
*/
	/*
	 	int arr_i;
	int abi;
	int ax; //y
	int bx; //x
	for(abi = 0; abi <= 1; abi++){
		//printf("BITMAP %d  {", abi);
	for(arr_i = 0;arr_i <= 2; arr_i++){
		//printf("ARRAY %d  (", arr_i);
	for(ax = 0; ax <= lim-1;ax++){
		for(bx = 0;bx <= lim-1;bx++){
			//if(image_test[abi][arr_i][ax][bx] != 0)
				printf("%d,", image_test[abi][arr_i][ax][bx]);
			//puts(image_test[abi][arr_i][ax][bx]);
		}
	}
	//printf(");");
	printf("|");
	}
	//printf("};\n");
	printf("#");
	}
	 */

	////
	//IPP
	//system("pause");
	//printf("ARR Time elapsed: %fs\n", ((double)clock() - start_arr) / CLOCKS_PER_SEC);
	return EXIT_SUCCESS;
}

/*
 * STANDARD FUNCTION
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
	char str[1024];
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
			+(image_in[i-1][j+1] * c)
			+(image_in[i][j-1] * c)
			+(image_in[i][j] * c)
			+(image_in[i][j+1] * c)
			+(image_in[i+1][j-1] * c)
			+(image_in[i+1][j] * c)
			+(image_in[i+1][j] * c);
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
	int centre_ax, centre_ay;
	int centre_bx, centre_by;
	int a, b;
	//int am, aa;
	//int bm, ba;
	//int cnt_array;
	double d[cnt_a][cnt_b];
	double dr[cnt_a][cnt_b];
	double dl[cnt_a][cnt_b];

	double stddev[cnt_a];
	double stddev_r[cnt_a];
	double stddev_l[cnt_a];

	double avg[cnt_a];
	double avg_r[cnt_a];
	double avg_l[cnt_a];

	double cnt_data = 0;
	int cnt_label = 0;
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
		size_a[a] = calc_size(image_label_in_a, a+L_BASE, &centre_ax, &centre_ay);
		length_a[a] = calc_length(image_label_in_a, a+L_BASE);
		ratio_a[a] = 4*PI*size_a[a]/(length_a[a]*length_a[a]);

		m=sprintf(&buf[posi], "SECTION A%d=>\n", a);posi += m;
		m=sprintf(&buf[posi], "%6d;|\n", size_a[a]);posi += m;
		for(b = 0; b < cnt_b; b++)
		{
			//m=sprintf(&buf[posi], "SECTION B%d=>\n", b);posi += m;
			size_b[b] = calc_size(image_label_in_b, b+L_BASE, &centre_bx, &centre_by);
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
				centre_ax, centre_ay, centre_bx, centre_by, a, b);
			//printf("\nCURRENT DISTANCE VALUE=> %f FROM %dA, %dB\n", d[a][b], a,b);
			//l[a][b] = abs(size_b[b]-size_a[a]);
			//m=sprintf(&buf[posi], ">DISTANCE:%d; => %d;<\n",b, d[a][b]);posi += m;
		}

		//END Ax[DATA] => Ax $$ Bx
		//SORT Ax[DATA] array

		qsort(d[a], cnt_b, sizeof(double), sort);
		qsort(dr[a], cnt_b, sizeof(double), sort);
		qsort(dl[a], cnt_b, sizeof(double), sort);
/*
		int z;
		printf("CURRENT ARRAY VALUE=> {");
		for(z = 0; z < cnt_b; z++)
			printf("%f, ", d[a][z]);
		//printf("}\n");
*/
		stddev[a] = stat_stddev(d[a], cnt_b);
		stddev_r[a] = stat_stddev(dr[a], cnt_b);
		stddev_l[a] = stat_stddev(dl[a], cnt_b);

		avg[a] = stat_avg(d[a], cnt_b);
		avg_r[a] = stat_avg(dr[a], cnt_b);
		avg_l[a] = stat_avg(dl[a], cnt_b);
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
					if(
							/*
							(
									(((d[a][0] - d[a][1]) >= stddev[a]) &&
									((d[a][0] - d[a][1]) > 0))

									||
									(((dl[a][0] - dl[a][1]) >= stddev_l[a]) &&
									((dl[a][0] - dl[a][1]) > 0))

									||
									(((dr[a][0] - dr[a][1]) >= stddev_r[a]) &&
									((dr[a][0] - dr[a][1]) > 0))
							)

							&&
							*/
							/*
							(
									((d[a][0] - avg[a]) > stddev[a])
									||
									((dl[a][0] - avg_l[a]) > stddev_l[a])
									||
									((dr[a][0] - avg_r[a]) > stddev_r[a])
							)
							*/

							((d[a][0]-d[a][0-1]) <= threshold_ff)

							&&
							((dl[a][0]-dl[a][0-1]) <= threshold_ff)

							&&
							((dr[a][0]-dr[a][0-1]) <= threshold_ff)


					)
					{
						//printf("PEAK!\n");
						m=sprintf(&buf[posi], ">PEAK DATA => %f\n", (d[a][0]-d[a][cnt_b-1]));posi += m;


						/*
						double cluster_value = 0;

						double cluster_avg = 0;
						double cluster_avg_l = 0;
						double cluster_avg_r = 0;


						double total_c = (double)cnt_b;

						int cpst_c = 0;
						int cpst_c_l = 0;
						int cpst_c_r = 0;

						double cpst_l[cnt_b];
						double cpst_l_l[cnt_b];
						double cpst_l_r[cnt_b];

						int compensate = 0;
						double cpst = 0;
						for(i = 0; i < cnt_b-1; i++)
						{
							if((d[a][i]-d[i][i+1]) > stddev[a])
							{
								cpst_l[cpst_c] = i;
								cpst_c++;
								compensate++;
							}
							if((dl[a][i]-dl[i][i+1]) > stddev_l[a])
							{
								cpst_l_l[cpst_c_l] = i;
								cpst_c_l++;
								//compensate++;
							}
							if((dr[a][i]-dr[i][i+1]) > stddev_r[a])
							{
								cpst_l_r[cpst_c_r] = i;
								cpst_c_r++;
								//compensate++;
							}

						}


						cluster_avg = stat_avg(cpst_l, cpst_c);
						cluster_avg_l = stat_avg(cpst_l, cpst_c_l);
						cluster_avg_r = stat_avg(cpst_l, cpst_c_r);
						//for(i = 0; i < cpst_c; i++)
							//if(cpst_l[cpst_c] < ceil((cnt_b*threshold_ff)))
						//CALC STAT AVG
						cluster_value = (((total_c - cluster_avg)/total_c)+
								((total_c - cluster_avg_l)/total_c)+
								((total_c - cluster_avg_r)/total_c))/3;
						//DEFINE STAT VALUE


						cpst = (double)(compensate)/(double)cnt_b;
						//DEFINE MAIN CNT VALUE, HIGHER = CLOSER

						cnt_data = cnt_data + ((double)cluster_value);
						printf("Compensate CNT :%f\n", cluster_value);
						//printf("Cluster value:%f\n", cluster_value);



						*/

						cnt_data++;//GET PEAK DATA; MARK AS FOUND....
						cnt_label++;
						label[cnt_label] = (a+L_BASE);
						//printf("MAKING LABEL PROG  %d=> %d\n",cnt_data , label[cnt_data]);

						m=sprintf(&buf[posi], ">PEAK SECTION => %d\n", a);posi += m;

						//CHECK IF ONE LABEL THRESHOLD HAS LARGER VALUE THAN OTHERS IN RESULT(PEAK TARGET DATA):3
						//IF UNDETECTED(GARBAGE DATA), MAIN CNT DECREMENT; BEFORE END, CHECK CNT VALUE,
						//IF CNT = ORIGINAL CNT; THEN => MATRIX A = MATRIX B
						//IF CNT < ORIGINAL CNT *BY VALUE*; THEN => MATRIX A < MATRIX B BY *VALUE*
						//IF CNT = 0; THEN => MATRIX A != MATRIX B
					}

					else /*if(

							(
									(((d[a][0] - d[a][1]) < stddev[a]) &&
									((d[a][0] - d[a][1]) > 0))
									||
									(((dl[a][0] - dl[a][1]) < stddev_l[a]) &&
									((dl[a][0] - dl[a][1]) > 0))

									||
									(((dr[a][0] - dr[a][1]) < stddev_r[a]) &&
									((dr[a][0] - dr[a][1]) > 0))
							)

							) */
					{
						/*
						double cluster_value = 0;
						double cluster_avg = 0;
						double total_c = (double)cnt_b;
						int cpst_c = 0;
						double cpst_l[cnt_b];
						int compensate = 0;
						double cpst = 0;
						for(i= 0; i < cnt_b-1; i++)
						{
							if((d[a][i]-d[i][i+1]) > stddev[a])
							{

								cpst_l[cpst_c] = i;
								cpst_c++;
								//IDEAL !(i-(i-1))
								//BAD DATA !(i)
								compensate++;
							}
						}

						cluster_avg = stat_avg(cpst_l, cpst_c);
						cluster_value = (total_c - cluster_avg)/(double)cnt_b;
						cpst = (double)(cnt_b - compensate)/(double)cnt_b;
						cnt_data = cnt_data + ((double)cpst*(double)cluster_value);
						printf("Compensate:%f\n", cpst);
						printf("Cluster value:%f\n", cluster_value);
						//HIGHER CNT = BETTER RESULT
						//LOWER CNT = LOWER ACCURACY

						 */
					}

					//*cc++;

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

int sort(const void *x, const void *y)
{
  //return (*(double*)x - *(double*)y);
	double xx = *(double*)x, yy = *(double*)y;
	  if (xx > yy) return -1;
	  if (xx < yy) return 1;
	  return 0;
}
double stat_stddev(double stat[], int stat_count)
{
	int	i;
	double avg_c = 0;
	double stat_c;
	double stddev;
	//CALC AVG
	for(i = 0; i < stat_count; i++)
		avg_c = avg_c + stat[i];

	avg_c = avg_c/stat_count;
	//CALC STDDEV
	for(i = 0; i < stat_count-1; i++)
		stat_c = pow(stat[i]-avg_c, 2);

	stddev = sqrt((stat_c/stat_count));

	return stddev;
}
double stat_avg(double stat[], int stat_count)
{
	int	i;
	double avg_c = 0;
	//CALC AVG
	for(i = 0; i < stat_count; i++)
		avg_c = avg_c + stat[i];

	avg_c = avg_c/stat_count;


	return avg_c;
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
	//LOOP MATRIX CONTENT
	//int x,y;
	//unsigned char image_buf_b[Y_SIZE][X_SIZE];
	//unsigned char image_buf_a[Y_SIZE][X_SIZE];
/*
	for (x = 0; x < Y_SIZE; x++)
		for (y = 0; y < X_SIZE; y++)
		{
			image_buf_a[y][x] = 0;
			image_buf_b[y][x] = 0;
		}
*/

	int total = 0;
	//SHIFT MATRIX TO CENTRE
	//int shift_ax, shift_ay;
	//int shift_bx, shift_by;
	//int aa,bb;
	//MIDPOINT = [X_SIZE/2][Y_SIZE/2]
	//shift_ax = ((X_SIZE/2) - centre_ax);
	//shift_ay = ((Y_SIZE/2) - centre_ay);
	//shift_bx = ((X_SIZE/2) - centre_bx);
	//shift_by = ((Y_SIZE/2) - centre_by);
	//printf("\n=>SHIFT XY => !");
	//printf("!LB => %d A(X, Y) => A(%d, %d)", label_a, centre_ax, centre_ay);
	//printf("!LB => %d B(X, Y) => B(%d, %d)<=",label_b, centre_bx, centre_by);
	int ex = X_SIZE*2;
	int ey = Y_SIZE*2;

	unsigned char space_a[ey][ex];
	unsigned char space_b[ey][ex];
	shift_centre(image_label_a, space_a, label_a, centre_ax, centre_ay);
	shift_centre(image_label_b, space_b, label_b, centre_bx, centre_by);
	//CENTRE POINT SHIFT END
	//LOOP MATRIX CONTENT
	int	i, j;
	for (i = 0; i < ey; i++)
	{
		for (j = 0; j < ex; j++)
		{
			//AND OPERATOR
			if ((image_buf_a[i][j] == (label_a+L_BASE)) && (image_buf_b[i][j] == (label_b+L_BASE)))
			{
				//INCREMENTAL
				total++;
				//AND OPERATOR MATRIX
				//image_label_and[i][j] = HIGH;

			}
			//TODO ADD OTHER FUNCTIONS....
			//


			//
			//if(image_buf_a[i][j] == label_a)
				//aa++;
			//if(image_buf_b[i][j] == label_b)
				//bb++;
		}
	}
	//write_bmp_mono(image_label_and, "kirino_intersect.bmp");
	//write_bmp_mono(image_buf_a, "kirino_intersecta.bmp");
	//write_bmp_mono(image_buf_b, "kirino_intersectb.bmp");
	//printf("[]=>TOTAL AB => %d!",total);
	//printf("!SIZE A => %d <=",size_a);
	//printf("!SIZE B => %d <=\n",size_b);
	//aa = 0;bb = 0;
	//write_bmp_mono(image_label_and, "kirino_tmp.bmp");
	//LOOP MATRIX CONTENT END

	double rtn = (double)total/((double)size_a+(double)size_b-(double)total);
	//printf("!RTN => %f!<=\n",rtn);
	return rtn;
	//i/(a+b-i)/
}
void shift_centre(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE],
	int label, int x, int y)
{
	int i, j;

	int shift_x = ((X_SIZE) - x);
	int shift_y = ((Y_SIZE) - y);
	for (i = 0; i < Y_SIZE*2; i++)
		{
			for (j = 0; j < X_SIZE*2; j++)
			{
				if(image_in[i][j] == (L_BASE+label))
				{
					//if((i+shift_y < Y_SIZE) && (j+shift_x < X_SIZE))
					image_out[i+shift_y][j+shift_x] = image_in[i][j];
				}

				//MOVE BY SHIFT X & Y MATRIX
			}
		}
}
long features_moment(unsigned char image_in[Y_SIZE][X_SIZE], int label, int x_shift, int y_shift)
{
	int i, j;
	long moment = 0;

	if((x_shift != 0) && (y_shift != 0))
	i = i+x_shift; j = j+y_shift;//SET CENTRE

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

				}

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

/*
 * NOISE END
 */

