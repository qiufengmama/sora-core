#include "image.h"
/* Support lib image related functions
  �摜�f�[�^�쐬
   width :�摜�̉���
   height:�摜�̏c��
   depth :�P��f������̃r�b�g��(8 or 24)
 */
ImageData* createImage(int width,int height,int depth)
{
  ImageData *newimg;
  int byte_per_pixel;

  if(width<0 || height<0) return NULL;
  if(depth!=8 && depth!=24) return NULL;  // �P��f������̃r�b�g��(8,24�ȊO�̓G���[)
  
  newimg=malloc(sizeof(ImageData));
  if(newimg==NULL) return NULL;
  // 1��f�i�[����̂ɕK�v�ȃo�C�g�������߂�
  byte_per_pixel=depth/8;
  // �摜�f�[�^���i�[����̂ɕK�v�ȃ��������m��
  newimg->pixels=malloc(sizeof(BYTE)*byte_per_pixel*width*height);
  if(newimg->pixels==NULL) {
    free(newimg);
    return NULL;
  }
  // �e�v���p�e�B�l��ݒ�
  newimg->width=width;
  newimg->height=height;
  newimg->depth=depth;
  return newimg;
}

/* List1-3
   �摜�f�[�^�̔p��
 */
void disposeImage(ImageData *img)
{
  if(img->pixels!=NULL) free(img->pixels);
  free(img);
  return;
}

/* List1-4
  �摜�f�[�^��̉�f�l���擾
  x,y ��f�̍��W
  pix ��f�l���i�[����
 */
int getPixel(ImageData *img,int x,int y,Pixel *pix)
{
  int ret=1;
  int adr;  // ��f�̉摜��̈ʒu
  int dep,val;
  BYTE *pixels;

  if(img==NULL) return -1;
  if(img->pixels==NULL) return -1;
  // �摜�O�̍��W���w�肳�ꂽ�ꍇ�̏����i�ł��߂��摜��̉�f���Q�Ƃ���j
  if(x<0) {
    x=0;
    ret=0;
  }
  if(x >= img->width ) {
    x=img->width -1;
    ret=0;
  }
  if(y<0) {
    y=0;
    ret=0;
  }
  if(y >= img->height ) {
    y=img->height -1;
    ret=0;
  }
  dep=img->depth;
  adr=x + y*img->width;
  pixels=img->pixels;
  if(dep==8) {  // �O���[�X�P�[���̏ꍇ�́ARGB���ׂĂ̓����l���Z�b�g����
    val=pixels[adr];
    pix->r=val;
    pix->g=val;
    pix->b=val;
  }
  else if(dep==24) {
    pixels+=(adr*3);
    pix->r=(*pixels);
    pixels++;
    pix->g=(*pixels);
    pixels++;
    pix->b=(*pixels);
  }
  else {
    return -1;
  }
  return ret;
}

/*
  ��f�l�̕␳�i�͈͊O�̒l��͈͓��Ɏ��߂�j
*/
int correctValue(int val,int max)
{
  if(val<0) return 0;
  if(val>max) return max;
  return val;
}

/* List1-5
  �摜�f�[�^��̉�f�l��ύX����
  x,y ��f�̍��W
  pix �Z�b�g�����f�l
 */
int setPixel(ImageData *img,int x,int y,Pixel *pix)
{
  int adr;  // ��f�̉摜��̈ʒu
  int dep,val;
  BYTE *pixels;

  if(img==NULL) return -1;
  if(img->pixels==NULL) return -1;
  // �摜�O�̍��W���w�肳�ꂽ��Ȃɂ����Ȃ�
  if(x<0 || x >= img->width || y<0 || y >= img->height ) {
    return 0;
  }
  dep=img->depth;
  adr=x + y*img->width;
  pixels=img->pixels;
  if(dep==8) {
    pixels[adr]=correctValue(pix->r,PIXELMAX);
  }
  else if(dep==24) {
    pixels+=(adr*3);
    (*pixels)=correctValue(pix->r,PIXELMAX);
    pixels++;
    (*pixels)=correctValue(pix->g,PIXELMAX);
    pixels++;
    (*pixels)=correctValue(pix->b,PIXELMAX);
  }
  else {
    return -1;
  }
  return 1;
}
