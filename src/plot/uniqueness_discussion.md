# 全ての結果

### 一意性

![ambiguous_prob.pdf](img/ambiguous_prob.pdf)
![multiple_prob.pdf](img/multiple_prob.pdf)
![noexist_prob.pdf](img/noexist_prob.pdf)

基本的な傾向

- n_obsが増加すると、一意性は上昇する
- Vmaxが増加すると、一意性は減少している（逆もあり）
- FOVが増加すると、一意性が減少している（逆もあり）
- thetaが増加すると、一意性が減少している
- kが増加すると、一意性が減少している

- n_obsが増加すると、一致しないケースが増える（正解が含まれない場合が増加するため）
- kが増加すると、一致しないケースが減る（正解が含まれない場合が減るため）



**傾向が全く逆になるケースが発生する理由を解明したい**

# V_mag = 1.5 の場合

![ambiguous_prob.pdf](img/Vmag1.5/ambiguous_prob.pdf)
![multiple_prob.pdf](img/Vmag1.5/multiple_prob.pdf)
![noexist_prob.pdf](img/Vmag1.5/noexist_prob.pdf) 

- どうも、FOVが増加すると一意性が増加しているケース（逆の状況）で発生しているような気がする。
  つまり、Vmagが小さいときに、FOVが小さくてn_obsが大きいとmultipleが増加している。


# FOVとmultiple_probの相関が-0.5以下の場合

![ambiguous_prob.pdf](img/FOVcorr<-0.5/ambiguous_prob.pdf)
![multiple_prob.pdf](img/FOVcorr<-0.5/multiple_prob.pdf)
![noexist_prob.pdf](img/FOVcorr<-0.5/noexist_prob.pdf)
![hist.pdf](img/FOVcorr<-0.5/hist.pdf)

- やっぱり、FOVとmultiple_probの相関が負の場合は、n_obsとmultiple_probの相関が正になっている。

# FOVとmultiple_probの相関が-0.5以上の場合

![ambiguous_prob.pdf](img/FOVcorr>-0.5/ambiguous_prob.pdf)
![multiple_prob.pdf](img/FOVcorr>-0.5/multiple_prob.pdf)![noexist_prob.pdf](img/FOVcorr>-0.5/noexist_prob.pdf)
![hist.pdf](img/FOVcorr>-0.5/hist.pdf)


やっぱり、FOVとmultiple_probの相関が負ではない場合は、n_obsとmultiple_probの相関が負になっている。

# なぜFOVとmultiple_probの相関が負になるのか

Vmax1.5theta_FOV0.08726646259971647 
![multiple_prob.pdf](img/Vmax1.5theta_FOV0.08726646259971647/multiple_prob.pdf) 
Vmax1.5theta_FOV0.5235987755982988 
![multiple_prob.pdf](img/Vmax1.5theta_FOV0.5235987755982988/multiple_prob.pdf) 
Vmax1.5theta_FOV0.17453292519943295 
![multiple_prob.pdf](img/Vmax1.5theta_FOV0.17453292519943295/multiple_prob.pdf) 
Vmax1.5theta_FOV1.0471975511965976 
![multiple_prob.pdf](img/Vmax1.5theta_FOV1.0471975511965976/multiple_prob.pdf) 

Vmax2.5theta_FOV0.08726646259971647 
![multiple_prob.pdf](img/Vmax2.5theta_FOV0.08726646259971647/multiple_prob.pdf) 
Vmax2.5theta_FOV0.5235987755982988 
![multiple_prob.pdf](img/Vmax2.5theta_FOV0.5235987755982988/multiple_prob.pdf) 
Vmax2.5theta_FOV0.17453292519943295 
![multiple_prob.pdf](img/Vmax2.5theta_FOV0.17453292519943295/multiple_prob.pdf) 
Vmax2.5theta_FOV1.0471975511965976 
![multiple_prob.pdf](img/Vmax2.5theta_FOV1.0471975511965976/multiple_prob.pdf)



Vmax5.5theta_FOV0.08726646259971647 
![multiple_prob.pdf](img/Vmax5.5theta_FOV0.08726646259971647/multiple_prob.pdf) 
Vmax5.5theta_FOV0.5235987755982988 
![multiple_prob.pdf](img/Vmax5.5theta_FOV0.5235987755982988/multiple_prob.pdf) 
Vmax5.5theta_FOV0.17453292519943295 
![multiple_prob.pdf](img/Vmax5.5theta_FOV1.0471975511965976/multiple_prob.pdf)  
Vmax5.5theta_FOV1.0471975511965976 
![multiple_prob.pdf](img/Vmax5.5theta_FOV1.0471975511965976/multiple_prob.pdf) 

 

- 連星（星の間隔＜theta_img）を観測する確率が増えてしまっているからかもしれない
  - Vmaxが小さいかつtheta_imgが大きい箇所でしか発生していない
  - Vmaxが大きいケースではn_obsが増加する従って微妙に増えているのも同じ理由かも
  - FOVが増えれば相対的に連星を観測する確率が小さくなるからかも









