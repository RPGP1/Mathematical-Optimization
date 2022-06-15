数理最適化特論 レポート１回目
===

###### tags: `数理最適化特論`

$$
\newcommand{\d}{\mathrm{d}}
\newcommand{\T}{\mathrm{T}}
\newcommand{\N}{\mathbb{N}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\C}{\mathbb{C}}
\newcommand{\{}{\left\{}
\newcommand{\}}{\right\}}
\newcommand{\(}{\left(}
\newcommand{\)}{\right)}
\newcommand{\<}{\left<}
\newcommand{\>}{\right>}
\newcommand{\vec}[1]{\mathbf{#1}}
\newcommand{\vecb}[1]{\boldsymbol{#1}}
\newcommand{\function}[1]{\,\mathrm{#1}\,}
\newcommand{\div}{\function{div}}
\newcommand{\rot}{\function{rot}}
\newcommand{\grad}{\function{grad}}
\newcommand{\diag}{\function{diag}}
\newcommand{\rank}{\function{rank}}
\newcommand{\Res}{\function{Res}}
\newcommand{\lcm}{\function{lcm}}
\newcommand{\Ker}{\function{Ker}}
\newcommand{\sgn}{\function{sgn}}
\newcommand{\tr}{\function{tr}}
\newcommand{\argmin}[1]{\underset{ #1 }{\arg\!\min}}
\newcommand{\argmax}[1]{\underset{ #1 }{\arg\!\max}}
\newcommand{\laplace}[2][]{\mathscr{L}^{ #1 }\left[ #2 \right]}
\newcommand{\ket}[1]{\left| #1 \>}
\newcommand{\bra}[1]{\< #1 \right|}
\newcommand{\bracket}[3][]{\< #2 \mid \! #1 \! \mid #3 \>}
\newcommand{\abs}[1]{\left| #1 \right|}
\newcommand{\norm}[1]{\left\| #1 \right\|}
\newcommand{\subs}[1]{\left. #1 \right|}
\newcommand{\diff}[2][]{
    \frac{\d #1}{\d #2}
}
\newcommand{\pdiff}[2][]{
    \frac{\partial #1}{\partial #2}
}
\newcommand{\align}[1]{
    \begin{align*}
        #1
    \end{align*}
}
\newcommand{array}[2][c]{
    \begin{array}{#1}
        #2
    \end{array}
}
\newcommand{\matrix}[2][c]{
    \left[
        \array[#1]{#2}
    \right]
}
\newcommand{\rmatrix}[2][c]{
    \left(
        \array[#1]{#2}
    \right)
}
$$

>   情報理工学系研究科 創造情報学専攻 修士2年 48-216622 秀島宇音
>   2022年 6月 13日

The source codes are available via [git](https://github.com/RPGP1/Mathematical-Optimization).

<!--
f = (Aw - b)^T (Aw - b) + l w^\T w
f = w^T (A^T A + l I) w - 2 b^\T A w + b^\T b

df = 2 (A^T A + l I) w - 2 A^T b
   = 2 (A^T A w + l w - A^T b)
ddf = 2 (A^T A + l I)
-->

## Q1

$A^\T A$の最大固有値を$\lambda_A$とすると、$f(\vec w)$は$(2(\lambda_A + \lambda))$-smoothである。そこで、

$$
\vec w_{k+1} = \vec w_k - \frac{\nabla f(\vec w_k)}{2(\lambda_A + \lambda)}
$$

という更新則の最急降下法を実装した。なお、$m = 3, n = 10$とした。

<center>

![](https://i.imgur.com/iuFE5oD.png)

<small>図: $\lambda = 0, 1, 10$のそれぞれにおける、イテレーション回数$k$と、真の最小値との誤差(倍精度浮動小数点)</small>
</center>

その結果が上の図である。全ての$\lambda$で一次収束に見える挙動が確認された。このことは次のように説明される。

*   $\lambda > 0$で$f(\cdot)$は強凸なので、最急降下法は一次収束になる。
*   $\lambda = 0$でも、$\R^n$における$\Ker A$の補空間では$f(\cdot)$は狭義凸であり、$\Ker A$に含まれる$\vec w$の成分は目的関数の値を変えないので、最急降下法は一次収束になる。

## Q2

$\lambda = 1$とし、Q1で実装した固定ステップの最急降下法と、新たに実装したArmijo's rule($\alpha_0 = 1,$ $\tau = 0.5,$ $\xi = 10^{-3}$)を適用する最急降下法とを比較した。

<center>

![](https://i.imgur.com/0iBeprQ.png)
<small>図: 固定ステップとArmijo's ruleの、イテレーション回数 -- 誤差の関係の比較</small>
</center>

まず、イテレーション回数$k$に対する誤差の減り方を比較したのが上の図である。どちらも一次収束だが、Armijo's ruleの方が各イテレーションでの誤差の減少率が大きかった。

|   |収束までの時間 \[ns\]|
|:-:|-:|
|固定ステップ|10,779|
|Armijo's rule|10,870|

<small>表: 固定ステップとArmijo's ruleの、収束までにかかる時間の比較</small>

次に、収束までにかかる時間を比較した結果が上の表である。その具体的な設定は以下の通り。

*   実行環境は次の通り
    *   CPU: Intel(R) Core(TM) i7-10750H CPU @ 2.60GHz
        *   ISA: x86_64
        *   L1Dキャッシュ: 192 KiB
        *   L1Iキャッシュ: 192 KiB
        *   L2キャッシュ: 1.5 MiB
        *   L3キャッシュ: 12 MiB
    *   メモリ: 32 GiB
    *   OS: Linux 5.14.0
*   単一スレッド実行
*   3回のwarm-upの後、10^16^回の計測の平均値を採用した
*   高速化のため、固定ステップの実装では10イテレーションごとに収束判定をした
    *   分岐の存在しない固定回数のループ文はloop-vectorizationが働いて高速化するため
*   収束判定は、前回の判定タイミングから$f(\cdot)$の値がどれだけ減少したのかを計算し、1イテレーションあたりの減少幅が$10^{-14}$より小さくなったら収束と判定した

Armijo's ruleによって実行時間は短縮されなかった。この理由として、最急降下法の計算が軽いものであるのに対し、Armijo's ruleによってループ・条件分岐のある更新則となることのオーバーヘッドが大きかったのではないかと推測する。

## Q3

$\lambda = 1$とし、Q1で実装した固定ステップの最急降下法と、新たに実装したNesterovの加速勾配法($\tilde\alpha_k = \frac{1}{L} = \frac{1}{2(\lambda_A + \lambda)},$ $\tilde\beta_k = \frac{k}{k+3}$)とを比較した。

<center>

![](https://i.imgur.com/N2ktaGv.png)
<small>図: 固定ステップ最急降下法とNesterovの加速勾配法の、イテレーション回数 -- 誤差の関係の比較</small>
</center>

まず、イテレーション回数$k$に対する誤差の減り方を比較したのが上の図である。どちらも一次収束だが、Nesterovの加速勾配法の方が各イテレーションでの誤差の減少率が小さかった。

|   |収束までの時間 \[ns\]|
|:-:|-:|
|固定ステップ最急降下法|10,977|
|Nesterovの加速勾配法|18,980|

<small>表: 固定ステップ最急降下法とNesterovの加速勾配法の、収束までにかかる時間の比較</small>

次に、収束までにかかる時間を比較した結果が上の表である。その具体的な設定は以下の点以外はQ2と同様である。

*   高速化のため、Nesterovの加速勾配法の実装では10イテレーションごとに収束判定をした
    *   分岐の存在しない固定回数のループ文はloop-vectorizationが働いて高速化するため

イテレーション回数の増え方と同程度に、Nesterovの加速勾配法は固定ステップ最急降下法と比べて実行時間が長くなった。
