centreXBegin, centreYBegin = 36.0, 42.0  #!!!!!!!!!!!!!!!!!position of centre of central zone at 24hrs!!!!!!!!!!!!!!!!
lowCutOffValue, highCutOffValue = -0.1, 0.5
diceThreshold = 0.0
meanNormalationForDataCleaning = False

L1centre =  [(71.693308612866261, 65.730379631939684, 17.508725039461556), (64.948373471739941, 74.85982027630439, 18.715769487444906), (69.626255034964515, 78.356788694860356, 17.510039930266206), (72.69216567759959, 83.446068169016954, 15.057104526160211), (42.661683909599908, 61.457352956002921, 16.777869572861992), (40.003016824085812, 55.425888791295158, 14.028463782436575), (35.668968322428647, 41.602230726685384, 24.283826671574072), (35.471906478248982, 47.25050133247808, 22.075485222450141), (28.476652294640171, 45.239348446726012, 21.185185416406338), (47.82499463069076, 53.752030175956698, 30.595295065226569), (44.512189280228945, 51.79770045505223, 34.70105823473363), (46.715365747832898, 48.568206138144291, 29.368057769535838), (47.716700322244769, 58.223058538368512, 31.538803662695667), (46.066876401869251, 54.036792452390308, 28.725288202041924), (49.944146807752922, 53.078896840120649, 26.739324105452386), (53.207006148602623, 50.185357709844112, 28.844991449369576), (51.179603047263889, 54.201333850693885, 32.803420639498967), (40.967546853393237, 48.004672253856427, 25.429751666447093), (53.731479922385667, 56.857677681683128, 25.948596081463503), (49.139571985682096, 54.522248619699702, 34.661061719756233)] #tracked backwards

#L1 cells within 30 microns of L1centre
meanVolumes = [129.97497878450014, 131.48852647720881, 132.5776909996965, 137.74725764424912, 141.45206204610645, 131.39257408805912, 145.01905658473197, 151.83856001755302, 155.27612169892828, 149.67308501547507, 158.5722402837913, 163.4224975424938, 161.59873688506721, 153.8707458905547, 164.8648864190879, 172.53397609493561, 163.42089130679796, 158.00999708684486, 166.21864803782645, 168.83368575076304]
#tracked backwards

centre3DBegin = [71, 67, 85*0.26/1.19]  #identified at 24hrs
indicesTBegin =  ["L1 central zone at 24hrs", [307, 315, 316, 321, 329, 331, 337, 341, 342, 344, 346, 349, 350, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 365, 366, 370, 371, 374, 375, 376, 377, 378, 379, 380, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435]]


#Parameters in identify indices script: [centreX, centreY, radius] = [54, 54, 30], zThr = zRes * 60, celldiff = 2  
indicesTEndBack = ["L1 central zone",  [808, 859, 860, 861, 618, 864, 875, 656, 676, 682, 685, 687, 693, 698, 885, 710, 801, 876, 887, 888, 734, 735, 755, 756, 765, 771, 772, 812, 778, 786, 708, 745, 793, 794, 798, 901, 807, 809, 813, 814, 733, 816, 818, 820, 822, 823, 827, 828, 829, 833, 835, 838, 844, 845, 848, 849, 850, 852, 654, 855, 856, 857, 858, 862, 863, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 877, 878, 879, 880, 881, 882, 884, 886, 889, 890, 891, 892, 893, 894, 895, 896, 898, 899, 904, 905, 906, 907, 900, 883, 832, 815, 764, 821, 841, 842, 843, 902, 761, 847, 897, 817, 851, 903]]

centre3DEnd = [48, 56, 156*0.26/1.127] 

indicesTEndForward = ['L1 full lineages after cleaning', [405, 457, 536, 569, 423, 487, 578, 468, 588, 539, 566, 552, 652, 683, 733, 563, 590, 648, 658, 710, 435, 519, 543, 613, 498, 677, 776, 521, 597, 662, 489, 649, 618, 693, 687, 461, 660, 483, 635, 614, 661, 682, 745, 786, 656, 654, 735, 734, 609, 685, 676, 708, 610, 589, 709, 502, 639, 638, 707, 727, 805, 752, 804, 388, 466, 592, 674, 628, 664, 740, 691, 741, 699, 724, 798, 772, 822, 684, 688, 779, 820, 787, 698, 771, 778, 813, 832, 723, 785, 722, 784, 718, 796, 744, 810, 831, 765, 756, 808, 818, 732, 795, 803, 824, 834, 853, 626, 743, 797, 770, 825, 841, 868, 867, 882, 807, 827, 835, 855, 737, 811, 773, 833, 843, 861, 842, 870, 753, 801, 761, 755, 764, 817, 816, 819, 830, 845, 812, 847, 809, 829, 783, 837, 836, 846, 838, 850, 862, 878, 884, 794, 793, 823, 851, 821, 856, 844, 869, 876, 828, 849, 857, 871, 860, 894, 877, 858, 864, 880, 889, 854, 865, 863, 885, 866, 891, 879, 904, 881, 893, 872, 859, 875, 886, 901, 890, 898, 852, 874, 888, 892, 895, 883, 887, 899, 900, 905, 906, 902, 903, 907, 896, 897] ]


#radius = 60, zThr = 10 , celldiff = 2  
indicesL1Init = [ "initial central L1", [332, 333, 340, 345, 346, 348, 349, 352, 353, 354, 357, 358, 363, 366, 367, 372, 374, 376, 377, 381, 383, 386, 387, 388, 390, 392, 393, 396, 398, 400, 401, 404, 406, 408, 409, 412, 415, 416, 417, 419, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488]]


