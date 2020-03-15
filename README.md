# A-convolution-kernel-implemented-by-Vivado-HLS
This project implements a convolution kernel based on vivado HLS on `zcu104`. All the hyper parameters
are fxed value, that is, M=16, OR=56, OC=56, N=16, IR=56, IC=56, K=3, S=1, P=1, so as to inform the speed up of our kernel compares to pure soft couterpart. Here OFM stands for output feature map, IFM stands for input feature map, W stands for weight, N stands for the number of input feature maps, m and n stands for the index of input feature map and output feature map, S stands for stride, r and c stands for row and column and P stands for Padding.
## Hardware struture
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

At the beginning of program firrst three rows and filter weight will be load to line buffer and filter
buffer (line13-14), then we process input data row by row (line15). The write row ifm will write one row
to line buffer at each iteration except the last iteration. The conv write function will process the input
line buffer (convolution operation) and write it to output buffer. The third function conv read will write
the output buffer to PS (PS is reading data from PL that is why it is called conv read). By carefully
write the code this three function can run simutaneously (line 21-22).
