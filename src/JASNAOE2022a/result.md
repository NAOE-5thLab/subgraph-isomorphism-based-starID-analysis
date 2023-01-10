

平均時間 ![time_mean.pdf](img/all/time_mean.pdf) 

- 天体の数と観測誤差に大きく依存していること
- 1sec以下で抑えるのは容易、Vmaxを5.5より小さくすれば良い



決定率

![determined_prob.pdf](img/ksplit/determined_prob.pdf)

 k < 1.0![determined_prob.pdf](img/k<1.0/determined_prob.pdf) 

k > 1.0

![determined_prob.pdf](img/k>1.0/determined_prob.pdf) 

k < 1 and n_obs > 2

![determined_prob.pdf](img/k<1andn_obs>2/determined_prob.pdf) 

k > 1 or n_obs = 2 ![determined_prob.pdf](img/k>1orn_obs=2/determined_prob.pdf) 

- kでおおきく傾向がへんかしていること
  - kが小さい
    - n_obsに反比例して性能が悪くなる
    - 天体の数にはほぼ左右されない
    - 観測誤差の大きさにも左右されない
  - kが大きい
    - n_obsは大きいほうが決定率が上昇
    - Vmaxは小さいほうが決定率が上昇
    - FOVは小さいほうが決定率が上昇
    - 観測誤差は小さいほうが決定率が上昇



観測確率

 ![obs_prob.pdf](img/k>1orn_obs=2/obs_prob.pdf) 

トレードオフ

![determined_prob.pdf](img/k>1orn_obs=2/determined_prob_obs_prob.pdf)

- Vmaxが小さい、FOVが小さいと観測確率が減少する
  - トレードオフの関係も見えるが両方とも高いケースも存在する
  - 極端な値を使用すればいいわけでもない





 正解率![correct_prob_on_determined.pdf](img/all/correct_prob_on_determined.pdf) 







 