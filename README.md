# A-convolution-kernel-implemented-by-Vivado-HLS
This project implements a convolution kernel based on `vivado HLS` on `zcu104`. All the hyper parameters
are fxed value, that is, M=16, OR=56, OC=56, N=16, IR=56, IC=56, K=3, S=1, P=1, so as to inform the speed up of our kernel compares to pure soft couterpart. Here OFM stands for output feature map, IFM stands for input feature map, W stands for weight, N stands for the number of input feature maps, m and n stands for the index of input feature map and output feature map, S stands for stride, r and c stands for row and column and P stands for Padding.

## Highlights
* A line buffer is used for input feature map buffer with Ping-Pong buffer protocol to reused and parallel process data.
* A Ping-Pong buffer is used for output buffer to maximise parallelism of program.
* The data pack method is used to increase the bandwidth of the whole system.
* The final speedup of this algorithm is 117.77x.

## Hardware struture
![hardware architecture](/fig/arch.png)
We identify data reuse opportunities in the feature maps of neighboring tiles. To this end, we naturally
implement line buffers. There are multiple channels of input feature maps (M) as shown in Figure 1.
Each line of the line buffers stores the same rows across all the channels. Winograd/Elemente-Wise-
Multiplication PEs fetch data from line buffers. We initiate an array of PEs by parallelizing the processing of the multiple channels. Finally, we use Ping-Pong buffers to overlap the data transfer and computation. All the input data (e.g. input feature maps,
filters) are stored in the external memory initially. The input and output feature maps are transferred
to FPGAs via a FIFO.

The first dimention of input line buffer is manually partioned to assist SDx recognize it as a rotate
buffer, moreover, it will help SDx to implement it as a Ping-Pong buffer. Thus we can write one rows to
input line buffer as the same time process three other input buffer rows. We process the input feature
maps row by row, this is implemented by a for loop (line1-4, line 20-21). Two output feature map buffers
are used to implement Ping-Pang buffer, for instance, we will firrst write data to buffer1 as the same time
read data from buffer2 in the next turn we will firrst write data to buffer2 as the same tiem read data
from buffer1 (line 21-22).

## Core code segment
At the beginning of program firrst three rows and filter weight will be load to line buffer and filter
buffer (line13-14), then we process input data row by row (line15). The write row ifm will write one row
to line buffer at each iteration except the last iteration. The conv write function will process the input
line buffer (convolution operation) and write it to output buffer. The third function conv read will write
the output buffer to PS (PS is reading data from PL that is why it is called conv read). By carefully
write the code this three function can run simutaneously (line 21-22).

```C
FIXDATA ifm_buff0[N][IC+2*P];
FIXDATA ifm_buff1[N][IC+2*P];
FIXDATA ifm_buff2[N][IC+2*P];
FIXDATA ifm_buff3[N][IC+2*P];

FIXDATA filter_buff[M][N][ROW_G][ROW_G];

FIXDATA ofm_buff0[M][OC];
FIXDATA ofm_buff1[M][OC];
int cifm_counter = 0;
int cofm_counter = 0;

load_cifm_data(cifm, ifm_buff1, ifm_buff2, ifm_buff3, &cifm_counter;
load_filter_buffer(tran_wgt, filter_buff);
for (int row = 0; row < IR; row += w_m)            
{

	if (rotate_counter == 0)
	{
		write_row_ifm(cifm, ifm_buff0, &cifm_counter, 1);
		conv_write(filter_buff, ifm_buff1, ifm_buff2, ifm_buff3, ofm_buff0);
		conv_read(cofm, ofm_buff1, &cofm_counter, row!=0);
	
	} else if (rotate_counter == 1)
	{
		write_row_ifm(cifm, ifm_buff1, &cifm_counter, 1);
		conv_write(filter_buff, ifm_buff2, ifm_buff3, ifm_buff0, ofm_buff1);
		conv_read(cofm, ofm_buff0, &cofm_counter, 1);
	} else if (rotate_counter == 2)
	{
		write_row_ifm(cifm, ifm_buff2, &cifm_counter, 1);
		conv_write(filter_buff, ifm_buff3, ifm_buff0, ifm_buff1, ofm_buff0);
		conv_read(cofm, ofm_buff1, &cofm_counter, 1);
	} else if (rotate_counter == 3)
	{
		write_row_ifm(cifm, ifm_buff3, &cifm_counter, row!=IR-1);
		conv_write(filter_buff, ifm_buff0, ifm_buff1, ifm_buff2, ofm_buff1);
		conv_read(cofm, ofm_buff0, &cofm_counter, 1);
	}
	rotate_counter += 1;
	if (rotate_counter == 4)
	{
		rotate_counter = 0;
	}
}
conv_read(cofm, ofm_buff1, &cofm_counter, 1);
```

## Experiment result
We have implemented four methods. In method one we using 16 Element-wise-multipliction PEs to calculate convolution without Ping-Pong buffer. The method two still utilized 16 Element-wise-multipliction PEs but using Ping-Pong of input and output buffer to acheive best parallelism. In method two although we can write, process convolution, write back data at the same time but the time of wrting one row is longer than process convolution, in another word our system is limited by the bandwith of PS-PL communication. Thus in method 3 we using data pack method to improve the system bandwidth by using 16 PEs. In method four we using the same method as method three except there are 32 PEs, in this method we improve the system computation ability to match the high bandwidth of current system and this acheive the best performance. The following compares the performance and resource utilization of these four methods.

![exp result](/fig/result.png)

## Reference
[1]	L. Lu, Y. Liang, Q. Xiao and S. Yan, "Evaluating Fast Algorithms for Convolutional Neural Networks on FPGAs," 2017 IEEE 25th Annual International Symposium on Field-Programmable Custom Computing Machines (FCCM), Napa, CA, 2017, pp. 101-108.

[2]	Shen, Junzhong \& Huang, You \& Wang, Zelong \& Qiao, Yuran \& Wen, Mei \& Zhang, Chunyuan. (2018). Towards a Uniform Template-based Architecture for Accelerating 2D and 3D CNNs on FPGA. 97-106. 10.1145/3174243.3174257. 

[3]	Chen Zhang, Peng Li, Guangyu Sun, Yijin Guan, Bingjun Xiao, and Jason Cong. 2015. Optimizing FPGA-based Accelerator Design for Deep Convolutional Neural Networks. In Proceedings of the 2015 ACM/SIGDA International Symposium on Field-Programmable Gate Arrays (FPGA '15). ACM, New York, NY, USA, 161-170. DOI: https://doi.org/10.1145/2684746.2689060

[4]	Guo, Kaiyuan \& Zeng, Shulin \& Yu, Jincheng \& Wang, Yu \& Yang, Huazhong. (2017). A Survey of FPGA Based Neural Network Accelerator. 
