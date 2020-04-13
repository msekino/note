# 数値計算
### Java
べき乗は Math.pow() を用いる。
```
Math.pow(2,10)
//1024
```
### Julia
TeXの記法と同じ ^ を用いる。
```
> 2^10
1024
```

# 配列の初期化
## 0で初期化した配列を作る
### Java
```
double[] x = new double[5];
double[][] X = new double[3][5];
```
### Julia
```
x = zeros(5)
X = zeros(3,5)
```

## 乱数がセットされた配列を得る
### Java
[0, 1)の一様乱数で一次元配列を生成
```
double[] x = new double[5];
Random rand = new Random();

for(int i = 0; i < x.length; i++) {
  x[i] = rand.nextDouble();
}
```
### Julia
[0, 1)の一様乱数で一次元配列を生成
```
x = rand(5)
```

## リストに値を貯めて配列に変換
### Java
Javaでは配列は静的配列なので、Array Listに突っ込んでから配列に変換する。  
ArrayListにはDoubleで入れる必要があるので、toArray()で得られるのもDoubleの配列になってしまう。  
プリミティブの配列に変換するには、さらに org.apache.commons.lang3.ArrayUtils.toPrimitive() を用いる。
```
List<Double> valueList = new ArrayList<>();

Random rand = new Random();

for(int i = 0; i < 5; i++) {
  valueList.add(rand.nextGaussian());
}

Double[] x2 = valueList.toArray(new Double[0]);
double[] x = ArrayUtils.toPrimitive(x2);
```

### Julia
Juliaの配列は実装は静的配列で高速なアクセスが可能でありつつ、動的配列のように値を追加することもできる。
```
x = []
for i in 1:5
    push!(x, randn())
end
```

# 行列演算
## 転置行列
### Java
```
public static double[][] transpose(double[][] X) {
  int dimRow = X.length;
  int dimCol = X[0].length;

  double[][] transpose = new double[dimCol][dimRow];

  for(int i = 0; i < dimRow; i++) {
    for(int j = 0; j < dimCol; j++) {
      transpose[j][i] = X[i][j];
    }
  }

  return transpose;
}
```
もし引数にfloatの配列を用いたいときには関数をさらに定義する必要がある。

### Julia
Juliaでは転置は'（プライム）で求まるが、もし関数として実装する場合、次のようになる。  
Anyのサブタイプ全ての配列に対して関数を定義できる。
特定の型に対して異なる挙動をさせたい場合には、その型に対する関数を定義すれば優先される。
```
function transpose(X::Array{<:Any,2})
    dimrow = size(X,1)
    dimcol = size(X,2)

    Y = Array{Any,2}(undef, dimcol, dimrow)

    for i in 1:dimrow
        for j in 1:dimcol
            Y[j,i] = X[i,j]
        end
    end

    return Y
end
```
戻り値の配列の型はArray{Any,2}になるが、値の型は入力値の型と同じになる。

## 逆行列、行列式
### Java
Javaでは逆行列や行列式は自力で実装するか、Javaのライブラリを利用する必要がある（使いやすいライブラリを知らない）。  
実装は手元にあるものの、かなりの行数になるので省略。
実装があれば、呼び出しはJuliaと同じくメソッド一行ではある。

### Julia
**逆行列**
```
inv(X)
```
**行列式**
```
det(X)
```
**トレース**
```
tr(X)
```

## 行列の積
### Java
```
public static double[][] product(double[][] X1, double[][] X2) {
  int dimRow = X1.length;
  int numMid = X2.length;
  int dimCol = X2[0].length;

  double[][] prods = new double[dimRow][dimCol];

  for(int i = 0; i < dimRow; i++) {
    for(int j = 0; j < dimCol; j++) {
      double prod = 0.0;

      for(int k = 0; k < numMid; k++) {
        prod += X1[i][k] * X2[k][j];
      }

      prods[i][j] = (double)prod;
    }
  }

  return prods;
}
```

### Julia
```
X * Y
```

## 要素ごとの積
### Java
```
public static void elementwiseProduct(double[][] X1, double[] X2) {
  double[][] prod = new double[X1.length][X1[0].length]

  for(int i = 0; i < X1.length; i++ ) {
    for(int j = 0; j < X1[0].length; j++) {
      prod[i][j] = X1[i][j] * X2[i][j];
    }
  }

  return prod;
}
```

### Julia
.（ドット演算子）でbroadcastして効率的に計算してくれる。
```
X .* Y
```

## 全ての要素に関数を適用
### Java
```
public static void apply(Function<Double, Double> f, double[][] X) {
  for(int i = 0; i < X.length; i++) {
    for(int j = 0; j < X[0].length; j++) {
      X[i][j] = f.apply(X[i][j]);
    }
  }
}
```
というメソッドを用意しておいて
```
apply(x -> Math.exp(x), X);
```
など。

### Julia
.（ドット演算子）でbroadcastして効率的に計算してくれる。
```
exp.(X)   #全ての要素にexp()を適用
```

# 機械学習の目的関数
## loglossの計算
確率の推定値 p の実現値 y ∈ {0,1} に対するloglossは
```
- y log(p) - (1-y) log(1-p)
```
だが、二回 log を計算するのを避けて

- y = 1 のとき -log(p)
- y = 0 のとき -log(1-p)

と実装する。
または一行で log も一回で -log((1-y) + (2y-1) * p) とも書ける。

### Java
Javaでは愚直に書くのが素直？  
zip関数が無いのでインデックスで参照する必要がある。
```
public double logloss(double[] p, double[] y) {
  double logloss = 0.0;

  for(int i = 0; i < y.length; i++) {
    if(y[i] == 1.0) {
      logloss += -Math.log(p[i]);
    } else {
      logloss += -Math.log(1.0-p[i]);
    }
  }

  return logloss;
}
```
または三項演算子を用いて
```
public double logloss(double[] p, double[] y) {
  double logloss = 0.0;

  for(int i = 0; i < y.length; i++) {
    logloss += (y[i]==1.0) ? -Math.log(p[i]) : -Math.log(1.0-p[i]);
  }

  return logloss;
}
```
など。

### Julia
Juliaにはzip関数があるのでfor文はインデックス参照より可読性が高くなる。
```
function logloss(ps, ys)
  logloss = 0.0;

  for (p, y) in zip(ps, ys)
    logloss += (y==1) ? -log(p) : -log(1-p)
  end

  return logloss;
end
```
内包表記を用いて
```
logloss(ps, ys) = sum([(y==1) ? -log(p) : -log(1-p) for (p, y) in zip(ps, ys)])
```
とも書ける。
一行の式を用いて
```
logloss(p, y) = - sum(@. log((1-y) + (2y-1) * p))
```
とも書ける。
@. は以降の演算を要素ごとに（broadcastして効率的に）演算するオペレータ。
