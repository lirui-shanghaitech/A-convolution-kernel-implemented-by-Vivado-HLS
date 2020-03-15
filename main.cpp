#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <cstring>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "sds_lib.h"


#define N 16
#define M 16
#define IR 56
#define IC 56
#define OR 56
#define OC 56
#define K 3
#define S 1
#define P 1
#define w_m 1
#define ROW_G 3
#define w_n 3  //w_n = w_m + w_r - 1 = 4 + 3 -1 = 6

typedef float  FIXDATA;

typedef struct ifm_struct
{
  float a0;
  float a1;
  float a2;
  float a3;
  float a4;
  float a5;
  float a6;
  float a7;
  float a8;
  float a9;
  float a10;
  float a11;
  float a12;
  float a13;
  float a14;
  float a15;
 } IPACK;

 typedef struct filter_struct
{
    float f0;
    float f1;
    float f2;
    float f3;
    float f4;
    float f5;
    float f6;
    float f7;
    float f8;
    float f9;
	float f10;
	float f11;
	float f12;
	float f13;
	float f14;
	float f15;

 }  FPACK;

 typedef struct ofm_struct
{
    float b0;
    float b1;
    float b2;
    float b3;
    float b4;
    float b5;
    float b6;
    float b7;
    float b8;
    float b9;
    float b10;
    float b11;
    float b12;
    float b13;
    float b14;
    float b15;
 }  OPACK;

void convolution_sw(float *ifm, float *ofm, float *weight)
{
    for (int m = 0; m < M; m++)
    {
        for (int r = 0; r < OR; r++)
        {
            for (int c = 0; c < OC; c++)
            {
                float odata = 0;
                int ofm_index = m*OR*OC + r*OC + c;
                for (int n = 0; n < N; n++)
                {
                    for(int kr = 0; kr < K; kr++)
                    {
                        for (int kc = 0; kc < K; kc++)
                        {
                            float ret;
                            int ic = c*S - P + kc;
                            int ir = r*S - P + kr;
                            int ifm_index = n*IR*IC + ir*IC + ic;
                            int wgt_index = m*N*K*K + n*K*K + kr*K + kc;

							if( (ic<0) || (ir<0) || (ic>(IC-1)) || (ir>(IR-1)))
                                ret=0;
                            else
                                ret = ifm[ifm_index];
                            ret *= weight[wgt_index];
                            odata += ret;
                        }
                    }
                }
                ofm[ofm_index] = odata;
            }
        }
    }
}

void generate(float* ifm, float* wgt)
{
	for(int i = 0; i < N*IR*IC; i++)
	{
        ifm[i] = (float)i/3000.0;
	}
	for(int i = 0; i < N*M*K*K; i++)
	{
        wgt[i] = (float)i/3000.0;
	}
}

void check(float* ofm_sw, float* ofm_hw)
{
	int error = 0;
	for(int i = 0; i < N*OR*OC; i++)
	{
		if (((ofm_sw[i] - ofm_hw[i]) > 0.01) || ((ofm_sw[i] - ofm_hw[i]) < -0.01))
			error++;
	}
	if(error>0) printf("error count = %d\n", error);
	else printf("correct!\n", error);
}


void change_ifm(float* ifm, IPACK* ifm_pack)
{
    //This function change the ifm to the standard form:
    /** LEN: #of input channels * #of elements for each row, Height: IR
     *  ROW ROW ROW .... ROW
     *  ROW ROW ROW .... ROW
     **/
    //Notice this part was implemented on PS, actually, in reality we don't need this part.
    //Since the input feature map is the output of previous feature, which is already the standard form.

    float ifm_temp[N][OR+2*P][OC+2*P];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < OR+2*P; j++)
        {
            for (int k = 0; k < OC+2*P; k++)
            {
                ifm_temp[i][j][k] = 0;
            }
        }
    }
    for (int bc = P; bc < OC+P; bc++)
    {
        for (int br = P; br < OR+P; br++)
        {
            for (int bn = 0; bn < N; bn++)
            {
                int temp1 = bn * IR * IC + (br-P)*IC;
                int i_index = temp1 + (bc-P);
                ifm_temp[bn][br][bc] = ifm[i_index];
            }
        }
    }

    int count = 0;
    for (int r = 0; r < OR+2*P; r++)
    {
        for (int c = 0; c < OR+2*P; c++)
        {

            ifm_pack[count].a0 = ifm_temp[0][r][c];
            ifm_pack[count].a1 = ifm_temp[1][r][c];
            ifm_pack[count].a2 = ifm_temp[2][r][c];
            ifm_pack[count].a3 = ifm_temp[3][r][c];
            ifm_pack[count].a4 = ifm_temp[4][r][c];
            ifm_pack[count].a5 = ifm_temp[5][r][c];
            ifm_pack[count].a6 = ifm_temp[6][r][c];
            ifm_pack[count].a7 = ifm_temp[7][r][c];
            ifm_pack[count].a8 = ifm_temp[8][r][c];
            ifm_pack[count].a9 = ifm_temp[9][r][c];
            ifm_pack[count].a10 = ifm_temp[10][r][c];
            ifm_pack[count].a11 = ifm_temp[11][r][c];
            ifm_pack[count].a12 = ifm_temp[12][r][c];
            ifm_pack[count].a13 = ifm_temp[13][r][c];
            ifm_pack[count].a14 = ifm_temp[14][r][c];
            ifm_pack[count].a15 = ifm_temp[15][r][c];
            count++;
        }
    }
}

void change_filter(float* wgt, FPACK* filter_pack)
{
    int count = 0;
    FIXDATA filter_buff[M][N][ROW_G][ROW_G];
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            for (int r = 0; r < ROW_G; r++)
            {
                for (int c = 0; c < ROW_G; c++)
                {
                    filter_buff[m][n][r][c] = wgt[count];
                    count++;
                }
            }
        }
    }
    count = 0;
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {

            filter_pack[count].f0 = filter_buff[m][n][0][0];
            filter_pack[count].f1 = filter_buff[m][n][0][1];
            filter_pack[count].f2 = filter_buff[m][n][0][2];
            filter_pack[count].f3 = filter_buff[m][n][1][0];
            filter_pack[count].f4 = filter_buff[m][n][1][1];
            filter_pack[count].f5 = filter_buff[m][n][1][2];
            filter_pack[count].f6 = filter_buff[m][n][2][0];
            filter_pack[count].f7 = filter_buff[m][n][2][1];
            filter_pack[count].f8 = filter_buff[m][n][2][2];
            count++;
        }
    }

}

void load_cifm_data(IPACK* cifm, FIXDATA ifm_buff0[N][(IC+2*P)], FIXDATA ifm_buff1[N][(IC+2*P)], FIXDATA ifm_buff2[N][(IC+2*P)], int* cifm_counter)
{
    /**
     * Load the standard input feature maps from PS to PL using FIFO, implemented at PL.
     * The ifm_buff is a cycle buffer, we put padding dirctly to the buffer. cifm_counter will
     * trace the current position of loaded feature map.
     **/
    for (int j = 0; j < (IC+2*P); j++)
    {
#pragma HLS PIPELINE
        ifm_buff0[0][j] = cifm[*cifm_counter].a0;
        ifm_buff0[1][j] = cifm[*cifm_counter].a1;
        ifm_buff0[2][j] = cifm[*cifm_counter].a2;
        ifm_buff0[3][j] = cifm[*cifm_counter].a3;
        ifm_buff0[4][j] = cifm[*cifm_counter].a4;
        ifm_buff0[5][j] = cifm[*cifm_counter].a5;
        ifm_buff0[6][j] = cifm[*cifm_counter].a6;
        ifm_buff0[7][j] = cifm[*cifm_counter].a7;
        ifm_buff0[8][j] = cifm[*cifm_counter].a8;
        ifm_buff0[9][j] = cifm[*cifm_counter].a9;
        ifm_buff0[10][j] = cifm[*cifm_counter].a10;
        ifm_buff0[11][j] = cifm[*cifm_counter].a11;
        ifm_buff0[12][j] = cifm[*cifm_counter].a12;
        ifm_buff0[13][j] = cifm[*cifm_counter].a13;
        ifm_buff0[14][j] = cifm[*cifm_counter].a14;
        ifm_buff0[15][j] = cifm[*cifm_counter].a15;
        (*cifm_counter)++;
    }


    for (int j = 0; j < (IC+2*P); j++)
    {
#pragma HLS PIPELINE
        ifm_buff1[0][j] = cifm[*cifm_counter].a0;
        ifm_buff1[1][j] = cifm[*cifm_counter].a1;
        ifm_buff1[2][j] = cifm[*cifm_counter].a2;
        ifm_buff1[3][j] = cifm[*cifm_counter].a3;
        ifm_buff1[4][j] = cifm[*cifm_counter].a4;
        ifm_buff1[5][j] = cifm[*cifm_counter].a5;
        ifm_buff1[6][j] = cifm[*cifm_counter].a6;
        ifm_buff1[7][j] = cifm[*cifm_counter].a7;
        ifm_buff1[8][j] = cifm[*cifm_counter].a8;
        ifm_buff1[9][j] = cifm[*cifm_counter].a9;
        ifm_buff1[10][j] = cifm[*cifm_counter].a10;
        ifm_buff1[11][j] = cifm[*cifm_counter].a11;
        ifm_buff1[12][j] = cifm[*cifm_counter].a12;
        ifm_buff1[13][j] = cifm[*cifm_counter].a13;
        ifm_buff1[14][j] = cifm[*cifm_counter].a14;
        ifm_buff1[15][j] = cifm[*cifm_counter].a15;
        (*cifm_counter)++;
    }

    for (int j = 0; j < (IC+2*P); j++)
    {
#pragma HLS PIPELINE
        ifm_buff2[0][j] = cifm[*cifm_counter].a0;
        ifm_buff2[1][j] = cifm[*cifm_counter].a1;
        ifm_buff2[2][j] = cifm[*cifm_counter].a2;
        ifm_buff2[3][j] = cifm[*cifm_counter].a3;
        ifm_buff2[4][j] = cifm[*cifm_counter].a4;
        ifm_buff2[5][j] = cifm[*cifm_counter].a5;
        ifm_buff2[6][j] = cifm[*cifm_counter].a6;
        ifm_buff2[7][j] = cifm[*cifm_counter].a7;
        ifm_buff2[8][j] = cifm[*cifm_counter].a8;
        ifm_buff2[9][j] = cifm[*cifm_counter].a9;
        ifm_buff2[10][j] = cifm[*cifm_counter].a10;
        ifm_buff2[11][j] = cifm[*cifm_counter].a11;
        ifm_buff2[12][j] = cifm[*cifm_counter].a12;
        ifm_buff2[13][j] = cifm[*cifm_counter].a13;
        ifm_buff2[14][j] = cifm[*cifm_counter].a14;
        ifm_buff2[15][j] = cifm[*cifm_counter].a15;
        (*cifm_counter)++;
    }
}

void write_row_ifm(IPACK* cifm, FIXDATA ifm_buff0[N][(IC+2*P)], int* cifm_counter, bool enable)
{
    if (enable)
    {

        for (int j = 0; j < (IC+2*P); j++)
        {
#pragma HLS PIPELINE
            ifm_buff0[0][j] = cifm[*cifm_counter].a0;
            ifm_buff0[1][j] = cifm[*cifm_counter].a1;
            ifm_buff0[2][j] = cifm[*cifm_counter].a2;
            ifm_buff0[3][j] = cifm[*cifm_counter].a3;
            ifm_buff0[4][j] = cifm[*cifm_counter].a4;
            ifm_buff0[5][j] = cifm[*cifm_counter].a5;
            ifm_buff0[6][j] = cifm[*cifm_counter].a6;
            ifm_buff0[7][j] = cifm[*cifm_counter].a7;
            ifm_buff0[8][j] = cifm[*cifm_counter].a8;
            ifm_buff0[9][j] = cifm[*cifm_counter].a9;
            ifm_buff0[10][j] = cifm[*cifm_counter].a10;
            ifm_buff0[11][j] = cifm[*cifm_counter].a11;
            ifm_buff0[12][j] = cifm[*cifm_counter].a12;
            ifm_buff0[13][j] = cifm[*cifm_counter].a13;
            ifm_buff0[14][j] = cifm[*cifm_counter].a14;
            ifm_buff0[15][j] = cifm[*cifm_counter].a15;
            (*cifm_counter)++;
        }
    }

}

void load_filter_buffer(FPACK* wgt, FIXDATA filter_buff[M][N][ROW_G][ROW_G])
{
    /**
     * Load the transformed filter weight to filter buffer. Implemented at PS.
     **/
    int count = 0;
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {

#pragma HLS PIPELINE
            filter_buff[m][n][0][0] = wgt[count].f0;
            filter_buff[m][n][0][1] = wgt[count].f1;
            filter_buff[m][n][0][2] = wgt[count].f2;
            filter_buff[m][n][1][0] = wgt[count].f3;
            filter_buff[m][n][1][1] = wgt[count].f4;
            filter_buff[m][n][1][2] = wgt[count].f5;
            filter_buff[m][n][2][0] = wgt[count].f6;
            filter_buff[m][n][2][1] = wgt[count].f7;
            filter_buff[m][n][2][2] = wgt[count].f8;
            count++;
        }
    }
}


void conv_write(FIXDATA filter_buff[M][N][ROW_G][ROW_G], FIXDATA ifm_buff0[N][(IC+2*P)], FIXDATA ifm_buff1[N][(IC+2*P)], FIXDATA ifm_buff2[N][(IC+2*P)], FIXDATA ofm_buff0[M][OC])
{
    for (int col = 0; col < IC; col += w_m)         //Column loop with stride w_m
    {
    #pragma HLS LOOP_TRIPCOUNT min=1 max=56
        for (int ti = 0; ti < M; ti++)              //Output channel loop
        {
	#pragma HLS PIPELINE
    #pragma HLS LOOP_TRIPCOUNT min=1 max=16
    #pragma HLS UNROLL FACTOR = 2
            FIXDATA Y = 0;
            for (int to = 0; to < N; to++)
            {
    #pragma HLS LOOP_TRIPCOUNT min=1 max=16
    #pragma HLS UNROLL FACTOR = 16
                FIXDATA mut000= ifm_buff0[to][col]*filter_buff[ti][to][0][0];
                FIXDATA mut100= ifm_buff1[to][col]*filter_buff[ti][to][1][0];
                FIXDATA mut200= ifm_buff2[to][col]*filter_buff[ti][to][2][0];
                FIXDATA mut010= ifm_buff0[to][col+1]*filter_buff[ti][to][0][1];
                FIXDATA mut110= ifm_buff1[to][col+1]*filter_buff[ti][to][1][1];
                FIXDATA mut210= ifm_buff2[to][col+1]*filter_buff[ti][to][2][1];
                FIXDATA mut020= ifm_buff0[to][col+2]*filter_buff[ti][to][0][2];
                FIXDATA mut120= ifm_buff1[to][col+2]*filter_buff[ti][to][1][2];
                FIXDATA mut220= ifm_buff2[to][col+2]*filter_buff[ti][to][2][2];
                FIXDATA acc000= mut000+mut100;
                FIXDATA acc010= mut200+mut010;
                FIXDATA acc020= mut110+mut210;
                FIXDATA acc030= mut020+mut120;
                FIXDATA acc040= acc000+acc010;
                FIXDATA acc050= acc020+acc030;
                FIXDATA acc060= acc040+acc050;
                Y += (acc060+mut220);
            }
            ofm_buff0[ti][col] = Y;
        }
    }
}

void conv_read(OPACK* cofm, FIXDATA ofm_buff0[M][OC], int* cofm_counter, bool enable)
{
    if (enable)
    {
        for (int j = 0;  j < OC; j++)
            {
#pragma HLS PIPELINE
                cofm[(*cofm_counter)].b0 = ofm_buff0[0][j];
                cofm[(*cofm_counter)].b1 = ofm_buff0[1][j];
                cofm[(*cofm_counter)].b2 = ofm_buff0[2][j];
                cofm[(*cofm_counter)].b3 = ofm_buff0[3][j];
                cofm[(*cofm_counter)].b4 = ofm_buff0[4][j];
                cofm[(*cofm_counter)].b5 = ofm_buff0[5][j];
                cofm[(*cofm_counter)].b6 = ofm_buff0[6][j];
                cofm[(*cofm_counter)].b7 = ofm_buff0[7][j];
                cofm[(*cofm_counter)].b8 = ofm_buff0[8][j];
                cofm[(*cofm_counter)].b9 = ofm_buff0[9][j];
                cofm[(*cofm_counter)].b10 = ofm_buff0[10][j];
                cofm[(*cofm_counter)].b11 = ofm_buff0[11][j];
                cofm[(*cofm_counter)].b12 = ofm_buff0[12][j];
                cofm[(*cofm_counter)].b13 = ofm_buff0[13][j];
                cofm[(*cofm_counter)].b14 = ofm_buff0[14][j];
                cofm[(*cofm_counter)].b15 = ofm_buff0[15][j];
                (*cofm_counter)++;
            }
    }

}

#pragma SDS data zero_copy(cifm[0:(IR+2*P)*(IC+2*P)], cofm[0:OR*OC], tran_wgt[0:N*M])
#pragma SDS data access_pattern(cifm: SEQUENTIAL, cofm: SEQUENTIAL, tran_wgt:SEQUENTIAL)
void convolution_hw(IPACK* cifm, OPACK* cofm, FPACK* tran_wgt)
{
    /**
     * In the first part, input feature maps, transformed weight maps will be loaded
     * from PS to PL using AXI interface (FIFO).
     **/
#pragma HLS data_pack variable=cifm struct_level
#pragma HLS data_pack variable=cofm struct_level
#pragma HLS data_pack variable=tran_wgt struct_level
    //Define ifm, ofm, wgt buffer and the partition way.
    FIXDATA ifm_buff0[N][IC+2*P];
#pragma HLS ARRAY_PARTITION variable=ifm_buff0 dim=1 complete
    FIXDATA ifm_buff1[N][IC+2*P];
#pragma HLS ARRAY_PARTITION variable=ifm_buff1 dim=1 complete
    FIXDATA ifm_buff2[N][IC+2*P];
#pragma HLS ARRAY_PARTITION variable=ifm_buff2 dim=1 complete
    FIXDATA ifm_buff3[N][IC+2*P];
#pragma HLS ARRAY_PARTITION variable=ifm_buff3 dim=1 complete

    FIXDATA filter_buff[M][N][ROW_G][ROW_G];
//#pragma HLS ARRAY_PARTITION variable=filter_buff dim=1 factor=4
#pragma HLS ARRAY_PARTITION variable=filter_buff dim=2 complete
#pragma HLS ARRAY_PARTITION variable=filter_buff dim=3 complete
#pragma HLS ARRAY_PARTITION variable=filter_buff dim=4 complete

    FIXDATA ofm_buff0[M][OC];
    FIXDATA ofm_buff1[M][OC];
#pragma HLS ARRAY_PARTITION variable=ofm_buff0 dim=1 factor=16
#pragma HLS ARRAY_PARTITION variable=ofm_buff1 dim=1 factor=16
    int cifm_counter = 0;
    int cofm_counter = 0;
    //These variable represent readable lines of line buffer.
    short unsigned int t0, t1, t2;
    //These varibale represent writable lines of line buffer.
    short unsigned int s1;
    short unsigned int rotate_counter = 0;

    /**Load data from PS. Here we can use dataflow pragma since filter weight and cifm can be load simutaneously.
     * load_filter_buffer may be one bottleneck of this design since we load all data from PS, actually, one can
     * load part of filter data by using bath method which beyond of this design's scope.
     **/

    load_cifm_data(cifm, ifm_buff1, ifm_buff2, ifm_buff3, &cifm_counter);
    load_filter_buffer(tran_wgt, filter_buff);

    /**
     * In the second part, caculating convolution using winograde method with line buffer for data reuse.
     **/
    //Define flag of ofm buffer, when flag = 0, writing to buffer1 and read from buffer2.
//#pragma HLS DATAFLOW
    for (int row = 0; row < IR; row += w_m)             //Row loop with stride w_m
    {
#pragma HLS LOOP_TRIPCOUNT min=1 max=16
        //In this subpart we need to define where is the legal line buffer, and rotate the line of line buffer.
        if (rotate_counter == 0)
        {
            write_row_ifm(cifm, ifm_buff0, &cifm_counter, 1);
            conv_write(filter_buff, ifm_buff1, ifm_buff2, ifm_buff3, ofm_buff0);
            conv_read(cofm, ofm_buff1, &cofm_counter, row!=0);

            // t0 = 1;
            // t1 = 2;
            // t2 = 3;
            // s1 = 0;
        } else if (rotate_counter == 1)
        {
            write_row_ifm(cifm, ifm_buff1, &cifm_counter, 1);
            conv_write(filter_buff, ifm_buff2, ifm_buff3, ifm_buff0, ofm_buff1);
            conv_read(cofm, ofm_buff0, &cofm_counter, 1);
            // t0 = 2;
            // t1 = 3;
            // t2 = 0;
            // s1 = 1;
        } else if (rotate_counter == 2)
        {
            write_row_ifm(cifm, ifm_buff2, &cifm_counter, 1);
            conv_write(filter_buff, ifm_buff3, ifm_buff0, ifm_buff1, ofm_buff0);
            conv_read(cofm, ofm_buff1, &cofm_counter, 1);
            // t0 = 3;
            // t1 = 0;
            // t2 = 1;
            // s1 = 2;
        } else if (rotate_counter == 3)
        {
            write_row_ifm(cifm, ifm_buff3, &cifm_counter, row!=IR-1);
            conv_write(filter_buff, ifm_buff0, ifm_buff1, ifm_buff2, ofm_buff1);
            conv_read(cofm, ofm_buff0, &cofm_counter, 1);
            // t0 = 0;
            // t1 = 1;
            // t2 = 2;
            // s1 = 3;
        }
        rotate_counter += 1;
        if (rotate_counter == 4)
        {
            rotate_counter = 0;
        }
    }
    conv_read(cofm, ofm_buff1, &cofm_counter, 1);
}

void chang_cofm(OPACK* ofm_pack, float* ofm)
{
    int count = 0;
    FIXDATA ofm_temp[M][OR][OC];
    int counter = 0;
    for (int r = 0; r < OR; r++)
    {
        for (int c = 0; c < OC; c++)
        {

                ofm_temp[0][r][c] = ofm_pack[counter].b0;
                ofm_temp[1][r][c] = ofm_pack[counter].b1;
                ofm_temp[2][r][c] = ofm_pack[counter].b2;
                ofm_temp[3][r][c] = ofm_pack[counter].b3;
                ofm_temp[4][r][c] = ofm_pack[counter].b4;
                ofm_temp[5][r][c] = ofm_pack[counter].b5;
                ofm_temp[6][r][c] = ofm_pack[counter].b6;
                ofm_temp[7][r][c] = ofm_pack[counter].b7;
                ofm_temp[8][r][c] = ofm_pack[counter].b8;
                ofm_temp[9][r][c] = ofm_pack[counter].b9;
                ofm_temp[10][r][c] = ofm_pack[counter].b10;
                ofm_temp[11][r][c] = ofm_pack[counter].b11;
                ofm_temp[12][r][c] = ofm_pack[counter].b12;
                ofm_temp[13][r][c] = ofm_pack[counter].b13;
                ofm_temp[14][r][c] = ofm_pack[counter].b14;
                ofm_temp[15][r][c] = ofm_pack[counter].b15;
                counter++;
        }
    }
    counter = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < OR; j++)
        {
            for (int k = 0; k < OC;  k++)
            {
                ofm[counter] = ofm_temp[i][j][k];
                counter++;
            }
        }
    }
}

void test_data_gen(float* ifm, float* wgt)
{
    printf("\nDebug: set ifm\n");
    int counter = 1;
    for (int bn = 0; bn < N; bn++)
    {
        for (int br = 0; br < OR; br++)
        {
            for (int bc = 0; bc < OC; bc++)
            {
                int temp1 = bn*IR * IC + br*IC;
                int i_index = temp1 + bc;
                ifm[i_index] = counter;
                counter++;
            }
        }
    }

    for (int i = 0; i < N*IR*IC; i++)
    {
        printf("%1f,", ifm[i]);
    }

    printf("\nDebug: set wgt\n");
    counter = 0;
    for (int wm = 0; wm < M; wm++)
    {
        for (int wn = 0; wn < N; wn++)
        {
            for (int wk1 = 0; wk1 < K; wk1++)
            {
                for (int wk2 = 0; wk2 < K; wk2++)
                {
                    int temp1 = wm*N * K * K + wn*K*K;
                    int temp2 = wk1*K + wk2;
                    int w_index = temp1 + temp2;
                    wgt[w_index] = counter;
                    counter++;
                }
            }
        }
    }
    for (int i = 0; i < N*M*K*K; i++)
    {
        printf("%1f,", wgt[i]);
    }
}

void printf_result(float* ofm, float* ofm_hw)
{
    printf("\nDebug: ofm result:\n");
    for (int i =  0; i < M*OR*OC; i++)
    {
        printf("[%1f, %1f]", ofm[i], ofm_hw[i]);
    }
}

int main()
{
	   float* ifm=(float*)sds_alloc(N*IR*IC*sizeof(float));
       FIXDATA* cifm=(FIXDATA*)sds_alloc(N*(IR+2*P)*(IC+2*P)*sizeof(FIXDATA));
	   float* ofm_sw=(float*)sds_alloc(M*OR*OC*sizeof(float));
	   float* ofm_hw=(float*)sds_alloc(M*OR*OC*sizeof(float));
       FIXDATA* cofm=(FIXDATA*)sds_alloc(M*OR*OC*sizeof(FIXDATA));
	   float* wgt   =(float*)sds_alloc(N*M*K*K*sizeof(float));
       IPACK *ifm_pack  = (IPACK *)(sds_alloc(sizeof(IPACK) * (IR+2*P)*(IC+2*P)));
       FPACK *filter_pack  = (FPACK *)(sds_alloc(sizeof(FPACK) * M*N));
       OPACK *ofm_pack  = (OPACK *)(sds_alloc(sizeof(OPACK) * (IR)*(IC)));
	   timeval start,end;

    //test_data_gen(ifm, wgt);

	printf("\nThis is EE216\n");

	generate(ifm,wgt);
	gettimeofday(&start, NULL);
	convolution_sw(ifm,ofm_sw,wgt);
	gettimeofday(&end, NULL);
	printf("\nconvolution_sw %lu us\n",(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));


    //First change ifm to standard input feature map ifm->cifm
    change_ifm(ifm, ifm_pack);
    change_filter(wgt, filter_pack);
    //Calculating the convolution and return cofm in standard form.
    gettimeofday(&start, NULL);
	convolution_hw(ifm_pack,ofm_pack,filter_pack);
	gettimeofday(&end, NULL);
	printf("\nconvolution_hw %lu us\n",(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));
    //Change the convolution result to required form. cofm->ofm_hw
    chang_cofm(ofm_pack, ofm_hw);
    //Finally, check if the data is right;
	check(ofm_sw,ofm_hw);
	return 0;
}
