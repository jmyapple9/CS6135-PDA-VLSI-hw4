# 最後實驗
<!-- ![image](https://hackmd.io/_uploads/H1AgD7F8p.png) -->
## epoch=3, numIter=100:50, lmbda=8000:8000, centerInit

testcase | wirelength |   runtime  | status
---------|------------|:----------:|:------:
public1  |  90456565  |      ?     | success
public2  |  12224783  |      ?     | success
public3  | 669056670  |      ?     | success
public7  | 244793337  |      ?     | success
public8  | 782512813  |      ?     | success
public9  | 809860885  |      ?     | success

## epoch=3, numIter=100:50, lmbda=4000:8000, centerInit

  testcase | wirelength |    runtime | status
-----------|------------|------------|--------
   public1 |   84470408 |      55.48 | success
   public2 |   16069494 |     123.45 | success
   public3 |  629132414 |     206.41 | success
   public7 |  246975589 |      78.75 | success
   public8 |  681815838 |     167.79 | success
   public9 |  832571331 |     210.02 | success

## epoch=3, numIter=100:50, lmbda=4000:4000, centerInit

testcase   | wirelength |    runtime | status
-----------|------------|:----------:|--------
public1    |   84112402 |      54.59 | success
public2    |   16792504 |     127.32 | success
public3    |  632482308 |     206.65 | success
public7    |  249702425 |      80.25 | success
public8    |  700911794 |     178.89 | success
public9    |  820234860 |     209.64 | success

## epoch=3, numIter=100:50, lmbda=4000:4000, randmoInit (left terminal)

exist NaN
testcase | wirelength |    runtime | status
-----------|------------|------------|--------
public1 |   94051135 |      51.77 | success
public2 |        N/A |        N/A | There is an error in the output results of public2 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
public3 |        N/A |        N/A | There is an error in the output results of public3 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
public7 |  320117506 |      78.41 | success
public8 |        N/A |        N/A | There is an error in the output results of public8 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
public9 |        N/A |        N/A | There is an error in the output results of public9 ([Error] Fail to Legalize! The global placement result cannot be legalized.).

## epoch=3, numIter=100:50, lmbda=2000:2000, centerInit (right terminal)

exist NaN
testcase | wirelength |    runtime | status
-----------|------------|------------|--------
public1 |   80872968 |      55.52 | success
public2 |   16044600 |     136.28 | success
public3 |  638797189 |     223.37 | success
public7 |        N/A |        N/A | There is an error in the output results of public7 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
public8 |        N/A |        N/A | There is an error in the output results of public8 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
public9 |        N/A |        N/A | There is an error in the output results of public9 ([Error] Fail to Legalize! The global placement result cannot be legalized.).


## epoch=3, numIter=100:50, lmbda=4000:2000, centerInit (best so far)

testcase   | wirelength |    runtime | status
-----------|------------|------------|--------
public1    |   83827612 |      56.46 | success
public2    |   16697557 |     127.17 | success
public3    |  633868222 |     207.51 | success
public7    |  249330319 |      80.23 | success
public8    |  700343528 |     177.15 | success
public9    |  824290117 |     212.94 | success

## epoch=3, numIter=100:50, lmbda=2000:4000, centerInit

exist Nan
testcase   | wirelength |    runtime | status
-----------|------------|------------|--------
public1    |   81145977 |      54.36 | success
public2    |   15631332 |     127.46 | success
public3    |  580706724 |     203.60 | success
public7    |        N/A |        N/A | There is an error in the output results of public7 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
public8    |        N/A |        N/A | There is an error in the output results of public8 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
public9    |        N/A |        N/A | There is an error in the output results of public9 ([Error] Fail to Legalize! The global placement result cannot be legalized.).


## epoch=3, numIter=100:50, lmbda=3000:2000, centerInit

exist Nan
  testcase | wirelength |    runtime | status
-----------|------------|------------|--------
   public1 |   79789534 |      51.09 | success
   public2 |   16708909 |     126.19 | success
   public3 |  540263794 |     211.78 | success
   public7 |        N/A |        N/A | There is an error in the output results of public7 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
   public8 |        N/A |        N/A | There is an error in the output results of public8 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
   public9 |  810440979 |     210.79 | success

## epoch=3, numIter=100:50, lmbda=4000:*1.5, centerInit

  testcase | wirelength |    runtime | status
-----------|------------|------------|--------
   public1 |   83784892 |      55.84 | success
   public2 |   16825777 |     134.19 | success
   public3 |  634253948 |     209.73 | success
   public7 |  248250911 |      82.67 | success
   public8 |  699602080 |     179.03 | success
   public9 |  825593831 |     204.06 | success


## epoch=3, numIter=100:50, lmbda=4000:0, stepsize=0.5, centerInit
  testcase | wirelength |    runtime | status
-----------|------------|------------|--------
   public1 |   84176662 |      50.72 | success
   public2 |   16811280 |     124.23 | success
   public3 |  642252834 |     207.74 | success
   public7 |        N/A |        N/A | There is an error in the output results of public7 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
   public8 |  724026217 |     176.31 | success
   public9 |  835291955 |     193.70 | success


## epoch=3, numIter=100:50, lmbda=4000:1000, centerInit
  testcase | wirelength |    runtime | status
-----------|------------|------------|--------
   public1 |   83828521 |      56.01 | success
   public2 |   16951158 |     136.23 | success
   public3 |  637555809 |     216.15 | success
   public7 |        N/A |        N/A | There is an error in the output results of public7 ([Error] Fail to Legalize! The global placement result cannot be legalized.).
   public8 |  711865925 |     186.36 | success
   public9 |  824826935 |     209.65 | success

## epoch=3, numIter=100:50, lmbda=4000:2000, centerInit, opt_gamma (left)
