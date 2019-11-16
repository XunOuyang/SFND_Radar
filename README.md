# SFND_Radar

**Implementation steps for the 2D CFAR process.**
The whole idea is similar to the 1D CFAR. 
![alt text](https://github.com/XunOuyang/SFND_Radar/blob/master/image/cfar2dbands.png)
*Image was from https://www.mathworks.com/help/phased/ref/2dcfardetector.html*
1. We save the bigger rectangle as z1. 
   The distance between the center point(target point) to the top of the rectange is Tr+Gr.
   The distance between the center point(target point) to the bottom of the rectange is Tr+Gr.
   So the total height of the rectangle is 1+2*(Tr+Gr).
   The distance between the center point(target point) to the left of the rectange is Td+Gd.
   The distance between the center point(target point) to the right of the rectange is Td+Gd.
   So the total height of the rectangle is 1+2*(Td+Gd).
2. We save the smaller rectangle as z2.
   The distance between the center point(target point) to the top of the rectange is Gr.
   The distance between the center point(target point) to the bottom of the rectange is Gr.
   So the total height of the rectangle is 1+2*Gr.
   The distance between the center point(target point) to the left of the rectange is Gd.
   The distance between the center point(target point) to the right of the rectange is Gd.
   So the total height of the rectangle is 1+2*Gd.
3. We sum the values within the bigger and smaller rectangle and convert the value from logarithmic to linear using db2pow function respetively. After the noise offset is added, we will get the final threshold. 
 

**Selection of Training, Guard cells and offset.**
1. From instinct, I select Tr=10, Td=8, Gr=6, Gd=3, offset=3.
    Based on rule of thumb, we do not usually use guard width 1, as it is too narrow for us. Especially the width of the whole data matrix is greater than 100. So 3 is selected here. As the guard width is 3, 1.5 times to 3 times of training cells will be an good option. As we have more rows than columns, so the 6 is taken as the Gr, and 12 is selected as Tr. offset 3 is selected from the previous quiz.
2. After we ran the code, we saw the result as below:
![alt text](https://github.com/XunOuyang/SFND_Radar/blob/master/image/1.PNG)
too much noise.

3. In order to remove the noise, we change the offset to 10.
![alt text](https://github.com/XunOuyang/SFND_Radar/blob/master/image/2.PNG)
we can get the range accurately. However the range of the range_rate is still too wide for us. We need to estimate the speed of the target accurately. So we have to adjust the Tr, Td, Gr, Gd a little bit more. Actually there could be more to be adjust for offset. As we increase the value of offset little by little, there doulbe multiple detections found. As the target shows on the figure would be clustered as 2 detections. Any values higher than 5 seems ok in this case.

4. When we increate the Td, Gd, the actual speed range increases. This is not what we want. 
Tr  Td  Gr  Gd  min_V max_V
12  8   6   3   17    25
12  12  6   6   15    26
20  8   8   4   15    25
We can take as many trials as we want. But Tr =12, Td = 8, Gr = 6, Gd = 3 seems to provide us a pretty decent performance.


**Steps taken to suppress the non-thresholded cells at the edges.**
All we need to do is just setting all the points which are out of range as 0. See the code as below:
```
RDM(union(1:(Tr+Gr),end-(Tr+Gr-1):end),:) = 0;  % Rows
RDM(:,union(1:(Td+Gd),end-(Td+Gd-1):end)) = 0;  % Columns 
```
