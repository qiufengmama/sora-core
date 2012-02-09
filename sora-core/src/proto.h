//DECL

//SORA CORE DECL

//STANDARD CONSTANT
int X_SIZE = 256;
int Y_SIZE = 256;
//DEBUG CONSTANT
char *debug_a = "lena.bmp";
char *debug_b = "lena_paint.bmp";
//lena_paint
//lena_mirror
//lena_colour


//STANDARD FUNCTION


char *join(const char* s1, const char* s2);
char *tochar(int var);
//STANDARD FUNCTION END

//STANDARD STAT FUNCTION

double stat_stddev(double stat[], int stat_count);
double stat_avg(double stat[], int stat_count);

//STANDARD STAT FUNCTION END

//SPECIAL STAT FUNCTION



//SPECIAL STAT FUNCTION END

//STANDARD DISTANCE DUNCTION

double euclidean_dist(double array[Y_SIZE*X_SIZE][2], int count);
//STANDARD DISTANCE FUNCTION END

void affine(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], 
	double a, double b, double c, double d, double e, double f);
void affine_coef(int x1, int y1, int u1, int v1, 
	int x2, int y2, int u2, int v2, 
	int x3, int y3, int u3, int v3,
	double *a, double *b, double *c, double *d, double *e, double *f);
void amplify(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double a, double b);
void back_subtraction(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_back[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int thresh);
void colorbar(unsigned char image_rgb[3][Y_SIZE][X_SIZE], int level);
void dilation(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE]);
void dither_ordered(unsigned char image_in[Y_SIZE][X_SIZE], 
		unsigned char image_out[Y_SIZE][X_SIZE]);
void dither_minimized(unsigned char image_in[Y_SIZE][X_SIZE], 
		unsigned char image_out[Y_SIZE][X_SIZE]);
void dither_minimized_multi(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int nq);
void dpcm1(unsigned char image_in[Y_SIZE][X_SIZE], int line,
	short data_out[X_SIZE]);
void dpcm2(unsigned char image_in[Y_SIZE][X_SIZE], int line,
	short data_out[X_SIZE]);
int dpcm_vlcode(unsigned char image_int[Y_SIZE][X_SIZE],
	unsigned char image_buft[Y_SIZE][X_SIZE]);
void erosion(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
int event(short dt);
void expand(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], int fmax, int fmin);
void expand_rgb(unsigned char image_in_rgb[3][Y_SIZE][X_SIZE], 
	unsigned char image_out_rgb[3][Y_SIZE][X_SIZE], int type);
void expand_ysh(int image_in_ysh[3][Y_SIZE][X_SIZE],
	int image_out_ysh[3][Y_SIZE][X_SIZE], int type);

/*
 * FEATURE
 */

/*
 * FFT
 */
int fft1(double a_rl[], double a_im[], int ex, int inv);
int fft2 (double a_rl[Y_SIZE][X_SIZE], double a_im[Y_SIZE][X_SIZE], int inv);
int fftfilter(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int a, int b);
int fftfilter_ff(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int a, int b);
int fftimage(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE]);
void fft1core(double a_rl[], double a_im[], int length, int ex, double sin_tbl[], double cos_tbl[], double buf[]);
void cstb(int length, int inv, double sin_tbl[], double cos_tbl[]);
void rvmtx1(double a[Y_SIZE][X_SIZE], double b[X_SIZE][Y_SIZE], int xsize, int ysize);
void rvmtx2(double a[X_SIZE][Y_SIZE], double b[Y_SIZE][X_SIZE], int xsize, int ysize);
void birv(double a[], int length, int ex, double b[]);
/*
 * FFT END
 */
/*
 * COLVOLUTION
 */
void convolute(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int b);
void assign(unsigned char image_out[Y_SIZE][X_SIZE], unsigned char image_in[Y_SIZE][X_SIZE]);
/*
 * CONVOLUTION END
 */
/*
 * THRESH
 */
void histgram(unsigned char image_in[Y_SIZE][X_SIZE], long hist[256]);
int threshdiscrim(long hist[256], double disparity);
void threshold(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int thresh, int type);
void thresh_color_difference(unsigned char image_in_rgb[3][Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int thresh, int type);
void threshold_discrim(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int type);
void threshold_dynamic(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int type);
unsigned char threshold_mode(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int smt, int type);
void thresh_rgb(unsigned char image_in_rgb[3][Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int thdrl, int thdrm, int thdgl, int thdgm, int thdbl, int thdbm);
void thresh_ysh(int image_in_ysh[3][Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE], int thdyl, int thdym, int thdsl, int thdsm, int thdhl, int thdhm);
/*
 * THRESH END
 */

/*
 * FEATURE
 */
int cmp(const void *x, const void *y);

double calc_distance(unsigned char image_label_a[Y_SIZE][X_SIZE],
	unsigned char image_label_b[Y_SIZE][X_SIZE],
	unsigned char image_label_and[Y_SIZE][X_SIZE],
	unsigned char image_buf_a[Y_SIZE][X_SIZE],
	unsigned char image_buf_b[Y_SIZE][X_SIZE],
	int size_a, int size_b,
	int centre_ax,  int centre_ay,
	int centre_bx,  int centre_by,
	int label_a, int label_b);

void features_compare(unsigned char	image_label_in_a[Y_SIZE][X_SIZE],
	unsigned char	image_label_in_b[Y_SIZE][X_SIZE],
	unsigned char image_label_out_a[Y_SIZE][X_SIZE],
	unsigned char image_label_out_b[Y_SIZE][X_SIZE],
	int cnt_a, int cnt_b, double size[], double length[],
	double ratio[],char *buf, double *rslt, int label[], int *cc, double threshold_ff);


double calc_size(unsigned char image_label[Y_SIZE][X_SIZE], int label, int *cx, int *cy);
double calc_equalize(unsigned char image_label[Y_SIZE][X_SIZE], int label, int *cx, int *cy);
double calc_length(unsigned char image_label[Y_SIZE][X_SIZE], int label);
double trace(unsigned char image_label[Y_SIZE][X_SIZE], int xs, int ys);
void extract_ratio(unsigned char image_label_in[Y_SIZE][X_SIZE],
	unsigned char image_label_out[Y_SIZE][X_SIZE],
	int cnt, double ratio[], double ratio_min, double ratio_max);
void extract_size(unsigned char image_label_in[Y_SIZE][X_SIZE],
	unsigned char image_label_out[Y_SIZE][X_SIZE],
	int cnt, double size[], double size_min, double size_max);
void features(unsigned char	image_label_in[Y_SIZE][X_SIZE],
	unsigned char image_label_out[Y_SIZE][X_SIZE],
	int cnt, double size[], double length[], double ratio[], char *buf);
// FEATURES COMPARE HERE...
double features_moment(unsigned char image_in[Y_SIZE][X_SIZE], int label);
void masking(unsigned char image_int[Y_SIZE][X_SIZE],
	unsigned char image_outt[Y_SIZE][X_SIZE],
	unsigned char image_maskt[Y_SIZE][X_SIZE]);
void label_masking(unsigned char image_lable[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE],int label[],int cc);
void labelset(unsigned char image[Y_SIZE][X_SIZE], int xs, int ys, int label);


double general_segmentation(unsigned char image_a[3][Y_SIZE][X_SIZE],unsigned char image_b[3][Y_SIZE][X_SIZE], int segment_size, int mode, double diffval[]);
double general_distance(unsigned char image_a[3][Y_SIZE][X_SIZE],unsigned char image_b[3][Y_SIZE][X_SIZE]);
void print_image(unsigned char image[3][Y_SIZE][X_SIZE],int channel);

/*
 * FEATURE END
 */

/*
 * NOISE
 */
void median(unsigned char image_in[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE]);
int median_value(unsigned char c[9]);

/*
 * NOISE END
 */


void frame_subtraction(unsigned char image_in1[Y_SIZE][X_SIZE],
	unsigned char image_in2[Y_SIZE][X_SIZE],
	unsigned char image_in3[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int thresh);
void gradient_difference(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double amp);
void gradient_roberts(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double amp);
void gradient_sobel(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double amp);
void hist2_image(unsigned char  image_in1[Y_SIZE][X_SIZE], 
	unsigned char image_in2[Y_SIZE][X_SIZE], 
	unsigned char image_hist[Y_SIZE][X_SIZE]);
void histimage(long hist[256], unsigned char image_hist[Y_SIZE][X_SIZE]);
void histprint(long hist[256], char *buf);
void histsmooth(long hist_in[256], long hist_out[256]);
void hough(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE],
	unsigned char image_hough[Y_SIZE][X_SIZE], int  thresh, char *buf);
void hough_cross(double rho1, double theta1, double rho2, double theta2,
	double *x, double *y);
void hue_image(int image_in_sat[Y_SIZE][X_SIZE],
	int image_in_hue[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int org);
void idpcm1(short data_in[X_SIZE], int line, 
	unsigned char image_out[Y_SIZE][X_SIZE]);
void idpcm2(short data_in[X_SIZE], int line, 
	unsigned char image_out[Y_SIZE][X_SIZE]);
int idpcm_vlcode(unsigned char image_buft[Y_SIZE][X_SIZE], 
	unsigned char image_outt[Y_SIZE][X_SIZE]);
int ievent(short ev);
void ivlcode(char vlc_in[], int no, short int data_out[]);
int labeling(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_label[Y_SIZE][X_SIZE]);
void laplacian(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double amp, int type);
void lattice(unsigned char image[Y_SIZE][X_SIZE]);
void log_zero_cross(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double var);


void mosaic(unsigned char image_in1[Y_SIZE][X_SIZE], 
	unsigned char image_in2[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE*2][X_SIZE*2],
	double mx, double my, double zm, double rt, int type);
void mosaic_affine(unsigned char image_in1[Y_SIZE][X_SIZE], 
	unsigned char image_in2[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE*2][X_SIZE*2],
	double a, double b, double c, double d, double e, double f, int type);
void mosaic_coef_blockmatch(unsigned char image_in1[Y_SIZE][X_SIZE], 
	unsigned char image_in2[Y_SIZE][X_SIZE],
	int x1, int y1, int x2, int y2, int xd, int yd, int bx, int by,
	int *mx, int *my);
void mosaic_coef_blockmatch_rgb(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	int x1, int y1, int x2, int y2, int xd, int yd, int bx, int by,
	int *mx, int *my);
void noise_rand(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int level);
void noise_spike(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int number, int level);
void perspect(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE],
	double a, double b, double c, double d, double e, double f,
	double g, double h);
void perspect_coef(int x1, int y1, int u1, int v1, 
	int x2, int y2, int u2, int v2, 
	int x3, int y3, int u3, int v3,
	int x4, int y4, int u4, int v4,
	double *a, double *b, double *c, double *d, double *e, double *f,
	double *g, double *h);
void plane(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], long hist[256]);
void prewitt(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double amp);
void pseudo_color(unsigned char image_in_m[Y_SIZE][X_SIZE],
	unsigned char image_out_r[Y_SIZE][X_SIZE], 
	unsigned char image_out_g[Y_SIZE][X_SIZE], 
	unsigned char image_out_b[Y_SIZE][X_SIZE], int type);
void quantize(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], int nq);
void radial_distortion(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double a, double b);
void range(unsigned char image_in[Y_SIZE][X_SIZE], int *fmax, int *fmin);
void rgb_to_ysh(unsigned char image_in_rgb[3][Y_SIZE][X_SIZE], 
	int image_out_ysh[3][Y_SIZE][X_SIZE]);
void rotation(unsigned char image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE], double deg);
void sat_image(int sat[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE]);
void scale(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double zx, double zy);
void scale_near(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double zx, double zy);
void scale_ng(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double zx, double zy);
void scale_rotate_shift(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], 
	double zx, double zy, double deg, double px, double py);
void shift(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], double px, double py);
void smooth(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], int type);
void smooth_edge_preserve(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE]);
void smooth_weighted(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE], int type);
void stereo_corre(unsigned char image_l[Y_SIZE][X_SIZE], 
	unsigned char image_r[Y_SIZE][X_SIZE], 
	unsigned char image_d[Y_SIZE][X_SIZE], int bsize);
void stereo_diff(unsigned char image_l[Y_SIZE][X_SIZE], 
	unsigned char image_r[Y_SIZE][X_SIZE], 
	unsigned char image_d[Y_SIZE][X_SIZE], int bsize);
void thinning(unsigned char	image_in[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
void tran_ysh(int image_in_ysh[3][Y_SIZE][X_SIZE],
	int image_out_ysh[3][Y_SIZE][X_SIZE], 
	double ya, double yb, double sa, double sb, double hb);
int vlcode(short int data_in[], int no, char vlc_out[]);
void y_image(int y[Y_SIZE][X_SIZE], unsigned char image_out[Y_SIZE][X_SIZE]);
void ysh_to_rgb(int image_in_ysh[3][Y_SIZE][X_SIZE],
	unsigned char image_out_rgb[3][Y_SIZE][X_SIZE]);
void zero_cross(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE]);
int check_bmp_file(char *filename);
int read_bmp_mono(unsigned char image[Y_SIZE][X_SIZE], char *filename);
int write_bmp_mono(unsigned char image[Y_SIZE][X_SIZE], char *filename);
int read_mono(unsigned char image[Y_SIZE][X_SIZE], char *filename);
int write_mono(unsigned char image[Y_SIZE][X_SIZE], char *filename);
int read_bmp_color(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int write_bmp_color(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int read_rgb_plane(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int write_rgb_plane(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int read_rgb_line(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int write_rgb_line(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int read_rgb_pixel(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int write_rgb_pixel(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int read_rgb_apart(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
int write_rgb_apart(unsigned char image[3][Y_SIZE][X_SIZE], char *filename);
void image_copy(unsigned char image_in[Y_SIZE][X_SIZE], 
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_copy_color(unsigned char image_in[3][Y_SIZE][X_SIZE], 
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
void image_clear(unsigned char image[Y_SIZE][X_SIZE]);
void image_clear_color(unsigned char image[3][Y_SIZE][X_SIZE]);
void image_negative(unsigned char image[Y_SIZE][X_SIZE]);
void image_negative_color(unsigned char image[3][Y_SIZE][X_SIZE]);
void color_to_mono(unsigned char image[3][Y_SIZE][X_SIZE]);
void image_addition(unsigned char image_in1[Y_SIZE][X_SIZE],
	unsigned char image_in2[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_addition_color(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
void image_subtraction(unsigned char image_in1[Y_SIZE][X_SIZE],
	unsigned char image_in2[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_subtraction_color(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
void image_difference(unsigned char image_in1[Y_SIZE][X_SIZE],
	unsigned char image_in2[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_difference_color(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
void image_multiplication(unsigned char image_in1[Y_SIZE][X_SIZE],
	unsigned char image_in2[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_multiplication_color(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
void image_or(unsigned char image_in1[Y_SIZE][X_SIZE],
	unsigned char image_in2[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_or_color(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
void image_and(unsigned char image_in1[Y_SIZE][X_SIZE],
	unsigned char image_in2[Y_SIZE][X_SIZE],
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_and_color(unsigned char image_in1[3][Y_SIZE][X_SIZE],
	unsigned char image_in2[3][Y_SIZE][X_SIZE],
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
void image_addition_rgb(unsigned char image[3][Y_SIZE][X_SIZE]);
void image_or_rgb(unsigned char image[3][Y_SIZE][X_SIZE]);
void image_and_rgb(unsigned char image[3][Y_SIZE][X_SIZE]);
void image_halfsize(unsigned char image_in[Y_SIZE*2][X_SIZE*2], 
	unsigned char image_out[Y_SIZE][X_SIZE]);
void image_halfsize_color(unsigned char image_in[3][Y_SIZE*2][X_SIZE*2], 
	unsigned char image_out[3][Y_SIZE][X_SIZE]);
 

