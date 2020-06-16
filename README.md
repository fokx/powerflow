# powerflow

A Power Flow Calculation Project,
 currently Newton-Raphson(N-R) Method and Fast Decoupled Power Flow(FDPF) Method implemented.
## Some examples
```shell script
$ python calculation_example1.py
```
### using example parameter set 1, visualized result
![Voltage distribution](network_voltage_distribution.png?raw=true "Title")

### using example parameter set 2 with N-R method

```shell script
Starting N-R calculation

delta P and delta Q array: before interation 1: 
[-2.0, 1.5, -1.0, 0.5]
max delta of abs deltaP & deltaQ: 2.0
***************Starting iteration 1**********************
----------
J
[[-20.  10.  -0.  -0.]
 [ 10. -20.  -0.  -0.]
 [  0.   0. -20.  10.]
 [  0.   0.  10. -20.]]
----------
dx
[[-0.08333333]
 [ 0.03333333]
 [-0.05      ]
 [-0.        ]]
----------
after num_iter: 1
V
[0.95, 1.0, 1.0]
theta(rad)
[-0.08333333333333333, 0.03333333333333334, 0.0]
theta(degree)
[-4.7746482927568605, 1.9098593171027445, 0.0]


delta P and delta Q array: before interation 2: 
[-0.10342852690588744, 0.06090762244129855, -0.14754650437400763, -0.07013451920063574]
max delta of abs deltaP & deltaQ: 0.14754650437400763
***************Starting iteration 2**********************
----------
J
[[-18.9024535    9.43542052   1.89657147   1.10582077]
 [  9.43542052 -19.42986548  -1.10582077  -1.43909238]
 [  1.89657147  -1.10582077 -17.1975465    9.43542052]
 [  1.10582077  -1.43909238   9.43542052 -20.57013452]]
----------
dx
[[-0.00693673]
 [ 0.00145902]
 [-0.01546035]
 [-0.0109761 ]]
----------
after num_iter: 2
V
[0.9353126693552949, 0.9890238998096266, 1.0]
theta(rad)
[-0.09027006707363609, 0.03479235623131796, 0.0]
theta(degree)
[-5.172093859682206, 1.9934551713702098, 0.0]


delta P and delta Q array: before interation 3: 
[-0.002966525718399282, 0.0020924019393127047, -0.0029324523650480216, -0.0008933784603666339]
max delta of abs deltaP & deltaQ: 0.002966525718399282
***************Starting iteration 3**********************
----------
J
[[-18.49326334   9.17821861   1.99703347   1.1538723 ]
 [  9.17821861 -19.06247211  -1.1538723   -1.4979076 ]
 [  1.99703347  -1.1538723  -16.49912824   9.17821861]
 [  1.1538723   -1.4979076    9.17821861 -20.06425887]]
----------
dx
[[-1.77943088e-04]
 [ 5.93248150e-05]
 [-3.17016627e-04]
 [-2.04204555e-04]]
----------
after num_iter: 3
V
[0.9350161596879176, 0.9888219366241157, 1.0]
theta(rad)
[-0.09044801016197766, 0.034851681046365095, 0.0]
theta(degree)
[-5.182289247637701, 1.9968542328928047, 0.0]


delta P and delta Q array: before interation 4: 
[-1.6315295905400262e-06, 1.2130534383647529e-06, -1.3795908273550594e-06, -3.5037233026002923e-07]
max delta of abs deltaP & deltaQ: 1.6315295905400262e-06

Calculation converges in required number of iterations(4 < 15).

Network loss: -4.1847615217527334e-07
Node 2 has maximum voltage of 1.0.
Node 0 has minimum voltage of 0.9350161596879176.
```
### using example parameter set 2 with FDPF
```shell script

Starting FDPF calculation

***************Starting iteration 1**********************
delta P array in interation 1: 
[0.145247175542351, -0.5, 0.11285171709832859]
delta P array in interation 1: 
[0.145247175542351, -0.5, 0.11285171709832859]
max delta P after abs: 0.5
max delta Q after abs: 0.5544570211602284
***************Starting iteration 2**********************
delta P array in interation 2: 
[0.06387362516582074, -0.055821219933934685, -0.025906860303969215]
delta P array in interation 2: 
[0.06387362516582074, -0.055821219933934685, -0.025906860303969215]
max delta P after abs: 0.06387362516582074
max delta Q after abs: 0.007221322832315746
***************Starting iteration 3**********************
delta P array in interation 3: 
[0.000744016403735337, -0.005743243995510072, -0.00296502569226359]
delta P array in interation 3: 
[0.000744016403735337, -0.005743243995510072, -0.00296502569226359]
max delta P after abs: 0.005743243995510072
max delta Q after abs: 0.0027685438611921853
***************Starting iteration 4**********************
delta P array in interation 4: 
[0.0009019499588351909, -0.0012420290695666636, -0.0012630719276042757]
delta P array in interation 4: 
[0.0009019499588351909, -0.0012420290695666636, -0.0012630719276042757]
max delta P after abs: 0.0012630719276042757
max delta Q after abs: 0.00022734393478973658
***************Starting iteration 5**********************
delta P array in interation 5: 
[0.0003404603899679337, -0.0002894434310528604, -0.00039591031615188177]
delta P array in interation 5: 
[0.0003404603899679337, -0.0002894434310528604, -0.00039591031615188177]
max delta P after abs: 0.00039591031615188177
max delta Q after abs: 6.118905989288548e-05
***************Starting iteration 6**********************
delta P array in interation 6: 
[7.595479483724077e-05, -5.915491373581494e-05, -0.000100128027527846]
delta P array in interation 6: 
[7.595479483724077e-05, -5.915491373581494e-05, -0.000100128027527846]
max delta P after abs: 0.000100128027527846
max delta Q after abs: 1.1626278852705507e-05
***************Starting iteration 7**********************
delta P array in interation 7: 
[1.4250356255463004e-05, -1.1257807754738725e-05, -2.3387310826017416e-05]
delta P array in interation 7: 
[1.4250356255463004e-05, -1.1257807754738725e-05, -2.3387310826017416e-05]
max delta P after abs: 2.3387310826017416e-05
max delta Q after abs: 2.0432969371486642e-06
***************Starting iteration 8**********************
delta P array in interation 8: 
[2.6124579193265163e-06, -2.12157502305077e-06, -5.340680202803316e-06]
delta P array in interation 8: 
[2.6124579193265163e-06, -2.12157502305077e-06, -5.340680202803316e-06]
max delta P after abs: 5.340680202803316e-06
max delta Q after abs: 3.759149135862394e-07

Calculation converges in required number of iterations(8 < 15) and epsilon 1e-05.

Network loss: 0.02768521615806724
Node 2 has maximum voltage of 1.05.
Node 0 has minimum voltage of 0.9695014886383434.
```

## comparison using example parameter set 1
```shell script
 $ time python calculation_example1.py 
```
N-R:

real    0m0.387s

user    0m0.714s

sys     0m2.225s

FDPF:

real    0m0.220s

user    0m0.380s

sys     0m1.177s

FDPF is faster than N-R method.