'''
# One line means one node record, blank lines does not matter

every line of nodes should of the order as below
[节点号，节点类型，发电机有功，发电机无功，负荷有功，负荷无功，节点电压初值
%      发电机有功下限   发电机有功上限    发电机无功下限   发电机无功上限    节点最低电压    节点最高电压]
'''

nodes2 = ''' 1   1   0.0000    0.0000  2.0000    1.0000        1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             2   1    1.5000    0.5000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             3   0    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
     '''
nodes1=''' 1   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             2   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             3   1    0.0000    0.0000    3.2200    0.0240    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             4   1    0.0000    0.0000    5.0000    1.8400    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             5   1       0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             6   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             7   1    0.0000    0.0000    2.3380    0.8400    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             8   1    0.0000    0.0000    5.2200    1.7600    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
             9   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            10   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            11   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            12   1    0.0000    0.0000    0.0750    0.8800    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            13   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            14   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            15   1    0.0000    0.0000    3.2000    1.5300    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            16   1    0.0000    0.0000    3.2940    0.3230    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            17   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            18   1    0.0000    0.0000    1.5800    0.3000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            19   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            20   1    0.0000    0.0000    6.8000    1.0300    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            21   1    0.0000    0.0000    2.7400    1.1500    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            22   1    0.0000    0.0000    0.0000    0.0000    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            23   1    0.0000    0.0000    2.4750    0.8460    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            24   1    0.0000    0.0000    3.0860   -0.9220    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            25   1    0.0000    0.0000    2.2400    0.4720    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            26   1    0.0000    0.0000    1.3900    0.1700    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            27   1    0.0000    0.0000    2.8100    0.7550    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            28   1    0.0000    0.0000    2.0600    0.2760    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            29   1    0.0000    0.0000    2.8350    0.2690    1.0000    0.00    0.00    0.00    0.00    0.95    1.07
            30  -1    2.5000    0.0000    0.0000    0.0000    1.0475    2.00    4.00    0.00    2.50    0.95    1.07
            31   0    0.0000    0.0000    0.0920    0.0460    0.9820    4.00    8.00    0.00    3.20    0.95    1.07
            32  -1    6.5000    0.0000    0.0000    0.0000    0.9831    4.00    8.00    0.00    3.20    0.95    1.07
            33  -1    6.3200    0.0000    0.0000    0.0000    0.9972    4.00    8.00    0.00    3.20    0.95    1.07
            34  -1    5.0800    0.0000    0.0000    0.0000    1.0123    4.00    8.00    0.00    3.20    0.95    1.07
            35  -1    6.5000    0.0000    0.0000    0.0000    1.0493    4.00    8.00    0.00    3.20    0.95    1.07
            36  -1    5.6000    0.0000    0.0000    0.0000    1.0635    4.00    8.00    0.00    3.20    0.95    1.07
            37  -1    5.4000    0.0000    0.0000    0.0000    1.0278    4.00    8.00    0.00    3.20    0.95    1.07
            38  -1    8.3000    0.0000    0.0000    0.0000    1.0265    5.00   10.00    0.00    3.60    0.95    1.07
            39  -1   10.0000    0.0000   11.0400    2.5000    1.0300    6.00   12.00    0.00    4.00    0.95    1.07
            '''


# i j R X half_Y
lines2='''
1 2 0 0.1 0 
1 3 0 0.1 0 
2 3 0 0.1 0 
'''
lines1='''
 1         2    0.0035   0.0411   0.34935
                 1        39    0.0010   0.0250   0.375
                 2         3    0.0013   0.0151   0.1286
                 2        25    0.0070   0.0086   0.073
                 3         4    0.0013   0.0213   0.1107
                 3        18    0.0011   0.0133   0.1069
                 4         5    0.0008   0.0128   0.0671
                 4        14    0.0008   0.0129   0.0691
                 5         6    0.0002   0.0026   0.0217
                 5         8    0.0008   0.0112   0.0738
                 6         7    0.0006   0.0092   0.0565
                 6        11    0.0007   0.0082   0.06945
                 7         8    0.0004   0.0046   0.039
                 8         9    0.0023   0.0363   0.1902
                 9        39    0.0010   0.0250   0.6
                10        11    0.0004   0.0043   0.03645
                10        13    0.0004   0.0043   0.03645
                13        14    0.0009   0.0101   0.08615
                14        15    0.0018   0.0217   0.183
                15        16    0.0009   0.0094   0.0855
                16        17    0.0007   0.0089   0.0671
                16        19    0.0016   0.0195   0.152
                16        21    0.0008   0.0135   0.1274
                16        24    0.0003   0.0059   0.034
                17        18    0.0007   0.0082   0.06595
                17        27    0.0013   0.0173   0.1608
                21        22    0.0008   0.0140   0.12825
                22        23    0.0006   0.0096   0.0923
                23        24    0.0022   0.0350   0.1805
                25        26    0.0032   0.0323   0.2565
                26        27    0.0014   0.0147   0.1198
                26        28    0.0043   0.0474   0.3901
                26        29    0.0057   0.0625   0.5145
                28        29    0.0014   0.0151   0.1245'''

transformers2='''
'''
# % Trans_para = [始端节点号 末端节点号  电阻   电抗  非标准变比]
transformers1='''
30         2     0.0000    0.0181    1.025
              31         6     0.0000    0.0250    1.07
              32        10     0.0000    0.0200    1.07
              11        12     0.0016    0.0435    1.006
              13        12     0.0016    0.0435    1.006
              33        19     0.0007    0.0142    1.07
              20        19     0.0007    0.0138    1.06
              34        20     0.0009    0.0180    1.009
              35        22     0.0000    0.0143    1.025
              36        23     0.0005    0.0272    1.0
              37        25     0.0006    0.0232    1.025
              38        29     0.0008    0.0156    1.025'''
parameter_set1 = (nodes1,lines1,transformers1)
parameter_set2 = (nodes2,lines2,transformers2)