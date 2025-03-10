// drawn wavetable from kahrs.pdf at 32k 200 Hz
 
const float table_kahrs000[160]  ={0.000000, 0.031403, 0.062805, 0.093933, 0.125275, 0.155975, 0.186798, 0.217163, 0.258942, 0.268646, 0.280518, 0.295013, 0.304535, 0.309357, 0.318939, 0.323730, 0.333344, 0.354919, 0.374115, 0.398071, 0.424438, 0.458069, 0.491577, 0.527588, 0.568359, 0.604279, 0.637939, 0.661835, 0.702637, 0.724243, 0.739075, 0.750610, 0.760742, 0.770081, 0.777802, 0.784698, 0.790131, 0.794434, 0.797577, 0.799347, 0.800018, 0.799347, 0.797607, 0.794403, 0.790192, 0.784576, 0.777924, 0.769958, 0.760864, 0.750549, 0.739075, 0.726562, 0.712738, 0.698090, 0.682007, 0.665283, 0.647095, 0.628357, 0.608246, 0.587494, 0.565704, 0.542969, 0.519653, 0.495148, 0.470398, 0.444305, 0.418121, 0.390778, 0.363281, 0.334869, 0.306213, 0.276855, 0.247192, 0.217224, 0.186615, 0.156250, 0.125000, 0.094147, 0.062683, 0.031433, 0.000000, -0.031372, -0.062805, -0.093994, -0.125183, -0.156067, -0.186737, -0.217163, -0.247223, -0.276855, -0.306213, -0.334839, -0.363312, -0.390778, -0.418091, -0.444397, -0.470276, -0.495239, -0.519562, -0.543030, -0.565735, -0.587402, -0.608398, -0.628143, -0.647552, -0.649872, -0.652252, -0.654724, -0.657013, -0.659515, -0.659454, -0.659454, -0.659546, -0.659363, -0.659576, -0.659393, -0.659515, -0.661865, -0.661865, -0.661865, -0.664307, -0.664215, -0.666687, -0.666687, -0.666626, -0.666748, -0.666595, -0.666656, -0.669128, -0.669006, -0.669098, -0.669067, -0.671448, -0.671448, -0.666718, -0.661804, -0.661926, -0.628235, -0.608307, -0.587494, -0.565674, -0.543030, -0.519592, -0.495239, -0.470245, -0.444458, -0.417999, -0.390869, -0.363251, -0.334869, -0.306183, -0.276886, -0.247192, -0.217163, -0.186768, -0.156036, -0.125183, -0.094025, -0.062775, 0.000031};

// drawn wavetable following fletcher: https://phys.unsw.edu.au/music/people/publications/Fletcher1988.pdf but at 44.1k and unsure of freq

const float table_fletcher000[99]  ={-0.853363, -0.853271, -0.853424, -0.866547, -0.906769, -0.906616, -0.906677, -0.906677, -0.906647, -0.893341, -0.893341, -0.880005, -0.879974, -0.880005, -0.880005, -0.880005, -0.880005, -0.880005, -0.879974, -0.880035, -0.879944, -0.880066, -0.879944, -0.880035, -0.880005, -0.866608, -0.866730, -0.866638, -0.866669, -0.853363, -0.853302, -0.853333, -0.853333, -0.839996, -0.840027, -0.839966, -0.826721, -0.826569, -0.813446, -0.799896, -0.800079, -0.786621, -0.786682, -0.773315, -0.746704, -0.733307, -0.693359, -0.613312, -0.493317, -0.333344, -0.160004, 0.013306, 0.173370, 0.306641, 0.413361, 0.519989, 0.586670, 0.653320, 0.733337, 0.773315, 0.800049, 0.813293, 0.840027, 0.813293, 0.826691, 0.813354, 0.813293, 0.813354, 0.840027, 0.826599, 0.813385, 0.799988, 0.786652, 0.773376, 0.746643, 0.719971, 0.680054, 0.639923, 0.586731, 0.519989, 0.439972, 0.373413, 0.333191, 0.226807, 0.173218, 0.106781, 0.053253, -0.053314, -0.133331, -0.186676, -0.266632, -0.320038, -0.386658, -0.466675, -0.519989, -0.586670, -0.626648, -0.653381, -0.693298};

// default from simforstacksansmem.c -> runsimplesir every 8 samples - somehow this needs normalising ! I tried...
const float plaguetable_simplesir[328]  ={-0.999996, -0.999996, -0.999996, -0.999996, -0.999995, -0.999995, -0.999994, -0.999994, -0.999993, -0.999993, -0.999992, -0.999992, -0.999991, -0.999991, -0.999990, -0.999989, -0.999988, -0.999987, -0.999986, -0.999985, -0.999984, -0.999983, -0.999981, -0.999980, -0.999978, -0.999977, -0.999975, -0.999973, -0.999971, -0.999968, -0.999966, -0.999963, -0.999960, -0.999957, -0.999954, -0.999950, -0.999946, -0.999942, -0.999938, -0.999933, -0.999927, -0.999922, -0.999916, -0.999909, -0.999902, -0.999894, -0.999886, -0.999877, -0.999867, -0.999857, -0.999846, -0.999834, -0.999821, -0.999807, -0.999792, -0.999775, -0.999758, -0.999739, -0.999718, -0.999696, -0.999672, -0.999647, -0.999619, -0.999589, -0.999557, -0.999522, -0.999485, -0.999445, -0.999401, -0.999354, -0.999304, -0.999249, -0.999191, -0.999127, -0.999059, -0.998985, -0.998906, -0.998820, -0.998728, -0.998628, -0.998521, -0.998405, -0.998281, -0.998146, -0.998001, -0.997845, -0.997676, -0.997494, -0.997298, -0.997087, -0.996859, -0.996614, -0.996349, -0.996064, -0.995756, -0.995424, -0.995067, -0.994682, -0.994266, -0.993818, -0.993336, -0.992815, -0.992255, -0.991650, -0.990999, -0.990297, -0.989541, -0.988726, -0.987847, -0.986901, -0.985881, -0.984783, -0.983599, -0.982325, -0.980952, -0.979473, -0.977880, -0.976165, -0.974319, -0.972331, -0.970191, -0.967888, -0.965410, -0.962743, -0.959874, -0.956788, -0.953470, -0.949902, -0.946066, -0.941944, -0.937515, -0.932758, -0.927650, -0.922167, -0.916283, -0.909971, -0.903203, -0.895949, -0.888177, -0.879855, -0.870948, -0.861420, -0.851235, -0.840354, -0.828737, -0.816344, -0.803132, -0.789060, -0.774085, -0.758165, -0.741255, -0.723315, -0.704303, -0.684179, -0.662906, -0.640449, -0.616775, -0.591858, -0.565674, -0.538205, -0.509439, -0.479372, -0.448008, -0.415357, -0.381440, -0.346288, -0.309943, -0.272455, -0.233888, -0.194316, -0.153824, -0.112508, -0.070475, -0.027842, 0.015266, 0.058714, 0.102363, 0.146066, 0.189675, 0.233035, 0.275997, 0.318407, 0.360118, 0.400986, 0.440874, 0.479652, 0.517198, 0.553402, 0.588163, 0.621394, 0.653016, 0.682967, 0.711194, 0.737658, 0.762331, 0.785196, 0.806248, 0.825491, 0.842939, 0.858615, 0.872549, 0.884776, 0.895340, 0.904289, 0.911674, 0.917550, 0.921976, 0.925010, 0.926715, 0.927151, 0.926382, 0.924469, 0.921473, 0.917456, 0.912476, 0.906592, 0.899859, 0.892333, 0.884066, 0.875109, 0.865509, 0.855315, 0.844571, 0.833318, 0.821598, 0.809448, 0.796905, 0.784003, 0.770776, 0.757253, 0.743465, 0.729438, 0.715197, 0.700767, 0.686171, 0.671430, 0.656564, 0.641591, 0.626529, 0.611395, 0.596203, 0.580968, 0.565703, 0.550421, 0.535133, 0.519849, 0.504581, 0.489337, 0.474126, 0.458955, 0.443833, 0.428765, 0.413759, 0.398821, 0.383955, 0.369166, 0.354460, 0.339840, 0.325311, 0.310875, 0.296536, 0.282297, 0.268160, 0.254129, 0.240204, 0.226389, 0.212684, 0.199092, 0.185613, 0.172250, 0.159002, 0.145871, 0.132858, 0.119963, 0.107187, 0.094530, 0.081993, 0.069575, 0.057277, 0.045099, 0.033041, 0.021102, 0.009283, -0.002417, -0.013999, -0.025461, -0.036806, -0.048033, -0.059143, -0.070136, -0.081014, -0.091776, -0.102423, -0.112957, -0.123377, -0.133685, -0.143880, -0.153965, -0.163940, -0.173805, -0.183561, -0.193210, -0.202752, -0.212188, -0.221519, -0.230745, -0.239868, -0.248889, -0.257807, -0.266626, -0.275344, -0.283964, -0.292486, -0.300910, -0.309239, -0.317472, -0.325611, -0.333657, -0.341611, -0.349472, -0.357244, -0.364925, -0.372518, -0.380023, -0.387441, -0.394773, -0.402019, -0.999997};

// also would be nice to have slower pulse - larger pause between pulses which we could maybe do in wavetable stuff 

// this one normalized in audacity but also speeded up so looks a bit rough...

const float plaguetable_simplesir_002[493]  ={-0.471802, -0.471741, -0.471619, -0.471558, -0.471466, -0.471313, -0.471283, -0.471130, -0.471039, -0.470917, -0.470795, -0.470673, -0.470581, -0.470398, -0.470276, -0.470154, -0.469971, -0.469849, -0.469696, -0.469513, -0.469360, -0.469147, -0.468994, -0.468842, -0.468536, -0.468475, -0.468109, -0.468018, -0.467743, -0.467438, -0.467316, -0.466919, -0.466797, -0.466461, -0.466125, -0.465942, -0.465515, -0.465332, -0.464844, -0.464600, -0.464172, -0.463898, -0.463440, -0.463165, -0.462585, -0.462341, -0.461761, -0.461395, -0.460876, -0.460419, -0.459869, -0.459381, -0.458801, -0.458191, -0.457764, -0.456940, -0.456451, -0.455750, -0.455048, -0.454437, -0.453583, -0.452881, -0.452087, -0.451294, -0.450378, -0.449646, -0.448486, -0.447876, -0.446564, -0.445831, -0.444611, -0.443604, -0.442413, -0.441345, -0.440063, -0.438812, -0.437622, -0.436066, -0.434906, -0.433289, -0.431915, -0.430267, -0.428772, -0.427002, -0.425385, -0.423553, -0.421753, -0.419769, -0.417877, -0.415771, -0.413666, -0.411530, -0.409180, -0.406952, -0.404388, -0.401978, -0.399353, -0.396698, -0.393860, -0.391144, -0.387909, -0.385132, -0.381775, -0.378571, -0.375214, -0.371643, -0.368103, -0.364288, -0.360504, -0.356506, -0.352417, -0.348145, -0.343811, -0.339203, -0.334686, -0.329803, -0.324829, -0.319824, -0.314423, -0.309113, -0.303467, -0.297729, -0.291840, -0.285614, -0.279541, -0.272827, -0.266449, -0.259430, -0.252472, -0.245270, -0.237793, -0.230225, -0.222412, -0.214325, -0.206207, -0.197662, -0.189087, -0.180115, -0.171143, -0.161713, -0.152374, -0.142487, -0.132538, -0.122467, -0.111908, -0.101440, -0.090515, -0.079407, -0.068207, -0.056641, -0.044922, -0.033081, -0.020844, -0.008575, 0.004059, 0.016663, 0.029755, 0.042694, 0.056213, 0.069489, 0.083344, 0.096954, 0.111084, 0.125122, 0.139374, 0.153870, 0.168213, 0.183075, 0.197540, 0.212585, 0.227386, 0.242340, 0.257507, 0.272430, 0.287659, 0.302795, 0.317902, 0.333130, 0.348236, 0.363373, 0.378448, 0.393494, 0.408478, 0.423370, 0.438171, 0.452850, 0.467529, 0.481903, 0.496368, 0.510406, 0.524567, 0.538330, 0.552032, 0.565552, 0.578735, 0.591919, 0.604614, 0.617310, 0.629700, 0.641754, 0.653748, 0.665131, 0.676727, 0.687622, 0.698456, 0.709045, 0.719025, 0.729248, 0.738617, 0.748077, 0.757080, 0.765717, 0.774231, 0.782257, 0.790100, 0.797577, 0.804749, 0.811676, 0.818115, 0.824524, 0.830292, 0.836151, 0.841309, 0.846527, 0.851196, 0.855682, 0.859924, 0.863739, 0.867462, 0.870758, 0.873840, 0.876648, 0.879120, 0.881531, 0.883453, 0.885345, 0.886719, 0.888214, 0.889099, 0.890106, 0.890594, 0.891022, 0.891266, 0.891113, 0.891052, 0.890350, 0.890106, 0.888916, 0.888214, 0.886902, 0.885559, 0.884125, 0.882416, 0.880554, 0.878662, 0.876404, 0.874298, 0.871674, 0.869385, 0.866425, 0.863892, 0.860748, 0.857819, 0.854645, 0.851318, 0.848053, 0.844513, 0.840942, 0.837341, 0.833466, 0.829773, 0.825684, 0.821808, 0.817657, 0.813538, 0.809296, 0.804993, 0.800598, 0.796234, 0.791718, 0.787201, 0.782623, 0.777924, 0.773285, 0.768524, 0.763733, 0.758942, 0.754028, 0.749146, 0.744263, 0.739197, 0.734283, 0.729218, 0.724121, 0.719147, 0.713898, 0.708923, 0.703613, 0.698547, 0.693298, 0.688171, 0.682831, 0.677704, 0.672302, 0.667236, 0.661896, 0.656586, 0.651428, 0.645966, 0.640839, 0.635437, 0.630127, 0.624908, 0.619568, 0.614227, 0.609009, 0.603577, 0.598389, 0.593048, 0.587799, 0.582458, 0.577209, 0.571899, 0.566681, 0.561371, 0.556183, 0.550781, 0.545715, 0.540405, 0.535217, 0.530029, 0.524811, 0.519623, 0.514496, 0.509277, 0.504150, 0.499084, 0.493866, 0.488831, 0.483734, 0.478638, 0.473572, 0.468567, 0.463501, 0.458496, 0.453552, 0.448486, 0.443604, 0.438629, 0.433685, 0.428772, 0.423920, 0.418945, 0.414246, 0.409271, 0.404510, 0.399719, 0.394897, 0.390167, 0.385406, 0.380646, 0.375977, 0.371307, 0.366577, 0.361969, 0.357391, 0.352631, 0.348236, 0.343506, 0.339081, 0.334503, 0.329926, 0.325531, 0.320953, 0.316650, 0.312073, 0.307739, 0.303345, 0.298920, 0.294678, 0.290253, 0.286011, 0.281677, 0.277405, 0.273224, 0.268890, 0.264832, 0.260468, 0.256439, 0.252197, 0.248199, 0.243958, 0.239960, 0.235870, 0.231812, 0.227875, 0.223724, 0.219910, 0.215729, 0.212036, 0.207886, 0.204163, 0.200104, 0.196442, 0.192383, 0.188782, 0.184784, 0.181091, 0.177338, 0.173553, 0.169861, 0.166168, 0.162384, 0.158844, 0.155060, 0.151581, 0.147888, 0.144287, 0.140747, 0.137177, 0.133667, 0.130096, 0.126678, 0.123047, 0.119873, 0.116058, 0.113037, 0.109222, 0.106293, 0.102448, 0.099579, 0.095764, 0.092957, 0.089172, 0.086395, 0.082611, 0.079926, 0.076141, 0.073456, 0.069824, 0.067017, 0.063507, 0.060760, 0.057281, 0.054474, 0.051208, 0.048218, 0.045197, 0.042023, 0.039276, 0.035889, 0.033447, 0.029816, 0.027649, 0.023773, 0.021973, 0.017792, 0.016388, 0.011841, 0.010895, 0.005890, 0.005493, 0.000061, 0.000031, -0.005554, -0.005432, -0.011108, -0.010864, -0.016571, -0.016266, -0.021820, -0.021820, -0.026825, -0.027588, -0.031433, -0.033630, -0.035431, -0.040527, -0.038147, -0.049164, -0.037598, -0.064240, -0.018890, -0.286316, -0.516144, -0.449127, -0.492249, -0.459412, -0.486328, -0.463409, -0.483368, -0.465546, -0.482880};

// this one is too fast at 48000 - create/test for 32000 - but what was source again - need better documentation
const float crowtable[142]  ={-0.283905, -0.055756, 0.206726, 0.495911, 0.636169, 0.630768, 0.330597, 0.458313, 0.337769, 0.321716, 0.154480, 0.151001, 0.140137, 0.057526, 0.132965, 0.294678, -0.039459, 0.053833, 0.190521, 0.194092, 0.362946, 0.368439, 0.407898, 0.434875, 0.309113, 0.305450, 0.217468, 0.098816, 0.233612, 0.127594, 0.140167, 0.386353, 0.420502, 0.355804, 0.382782, 0.267731, 0.079102, 0.194031, 0.197723, 0.165314, 0.167114, -0.055725, 0.073730, 0.129303, 0.195984, 0.100525, 0.095306, -0.145569, -0.215668, -0.043091, 0.001770, -0.073669, -0.307281, -0.402496, -0.672150, -0.557007, -0.722443, -0.857147, -0.867981, -0.882294, -0.880554, -0.803284, -0.724121, -0.632690, -0.438293, -0.330811, -0.224487, -0.201385, -0.258698, -0.249817, -0.246155, -0.226471, -0.422241, -0.217529, -0.096924, -0.093536, 0.079132, 0.141907, 0.165344, 0.019806, 0.159882, 0.079102, 0.267761, 0.170654, 0.203156, -0.097137, 0.070160, 0.273102, 0.061127, 0.265930, 0.328857, 0.546295, 0.864349, 0.828430, 0.654144, 0.569550, 0.469147, 0.510223, 0.413422, 0.567810, 0.539062, 0.776398, 0.600067, 0.573395, 0.468903, 0.316345, 0.098785, 0.077271, 0.077332, -0.050385, -0.312592, -0.479919, -0.334137, -0.016266, 0.131287, -0.143829, -0.366577, -0.637909, -0.513977, -0.497772, -0.368317, -0.334381, -0.499390, -0.469208, -0.303558, -0.323547, -0.023285, -0.160004, 0.079132, 0.053802, -0.044830, -0.172546, -0.161713, -0.123993, -0.346863, -0.244354, -0.447479, -0.587585, -0.460052, -0.598450, -0.321564, -0.025269};

// 96000 in audacity
const float crowtable_slower[283]  ={-0.264069, -0.214813, -0.075928, 0.079102, 0.227264, 0.363129, 0.475250, 0.570251, 0.656921, 0.690613, 0.610107, 0.452240, 0.351227, 0.375122, 0.437897, 0.427887, 0.357727, 0.313049, 0.302124, 0.259338, 0.173584, 0.116852, 0.132568, 0.169739, 0.157806, 0.095215, 0.040466, 0.054413, 0.148926, 0.262115, 0.279449, 0.149475, -0.025208, -0.074005, 0.040771, 0.176270, 0.202698, 0.162476, 0.183044, 0.285461, 0.372925, 0.382507, 0.359558, 0.371063, 0.415619, 0.445129, 0.428131, 0.373138, 0.314880, 0.290344, 0.300751, 0.293213, 0.221130, 0.123932, 0.095825, 0.163635, 0.235718, 0.215637, 0.126373, 0.078918, 0.140747, 0.272003, 0.386292, 0.433899, 0.419861, 0.381287, 0.356812, 0.364136, 0.381439, 0.359253, 0.269470, 0.147552, 0.077148, 0.108978, 0.196198, 0.236664, 0.195679, 0.147552, 0.167542, 0.208954, 0.164948, 0.035248, -0.053680, -0.019226, 0.071747, 0.119781, 0.130981, 0.162933, 0.194336, 0.165161, 0.101990, 0.084778, 0.094147, 0.020569, -0.144623, -0.257141, -0.216309, -0.103516, -0.042480, -0.028900, 0.001465, 0.013397, -0.073578, -0.220367, -0.307373, -0.328308, -0.402618, -0.561615, -0.671906, -0.637360, -0.557220, -0.586975, -0.722168, -0.832794, -0.857330, -0.852325, -0.867706, -0.885040, -0.882294, -0.878571, -0.880676, -0.858887, -0.803162, -0.749420, -0.724609, -0.698822, -0.632141, -0.531616, -0.439087, -0.376282, -0.329803, -0.279419, -0.225647, -0.192291, -0.199860, -0.236206, -0.260193, -0.255524, -0.248108, -0.255005, -0.248047, -0.216370, -0.224640, -0.322296, -0.424255, -0.385986, -0.215576, -0.086578, -0.098724, -0.145020, -0.091705, 0.020508, 0.077545, 0.090302, 0.143219, 0.204315, 0.164246, 0.051361, 0.020599, 0.105011, 0.159546, 0.109558, 0.078888, 0.168182, 0.268311, 0.245758, 0.169525, 0.175934, 0.204742, 0.095398, -0.099365, -0.131500, 0.072937, 0.279877, 0.269684, 0.123444, 0.065033, 0.151886, 0.261383, 0.308807, 0.333923, 0.404572, 0.540741, 0.718445, 0.870361, 0.909027, 0.822052, 0.709930, 0.660858, 0.636047, 0.562531, 0.480133, 0.476227, 0.520508, 0.502869, 0.429596, 0.420532, 0.503632, 0.560608, 0.538208, 0.545868, 0.661316, 0.769684, 0.733215, 0.606232, 0.545227, 0.567596, 0.559021, 0.473938, 0.381195, 0.311890, 0.218872, 0.102539, 0.044678, 0.074524, 0.108673, 0.079285, 0.014740, -0.051270, -0.154938, -0.312683, -0.447113, -0.478790, -0.425385, -0.336426, -0.200623, -0.012939, 0.132416, 0.126740, -0.002716, -0.138306, -0.241974, -0.373352, -0.535156, -0.630157, -0.601715, -0.522675, -0.487610, -0.488037, -0.456207, -0.379028, -0.314514, -0.323029, -0.408722, -0.511627, -0.539185, -0.456390, -0.345306, -0.317017, -0.350189, -0.309631, -0.157440, -0.037506, -0.067627, -0.145386, -0.094635, 0.064423, 0.143341, 0.068604, -0.033264, -0.059601, -0.076294, -0.157959, -0.225433, -0.176208, -0.082672, -0.109985, -0.260742, -0.360565, -0.310760, -0.231079, -0.286316, -0.460144, -0.588745, -0.575531, -0.493347, -0.471558, -0.535767, -0.587677, -0.519745, -0.331421, -0.129639, -0.015930};

// sine
const float ourtable[512]={0.000000, 0.012296, 0.024589, 0.036879, 0.049164, 0.061441, 0.073708, 0.085965, 0.098208, 0.110437, 0.122649, 0.134842, 0.147016, 0.159166, 0.171293, 0.183394, 0.195467, 0.207511, 0.219523, 0.231502, 0.243446, 0.255353, 0.267222, 0.279050, 0.290836, 0.302578, 0.314275, 0.325923, 0.337523, 0.349071, 0.360567, 0.372008, 0.383393, 0.394720, 0.405988, 0.417194, 0.428337, 0.439415, 0.450426, 0.461370, 0.472244, 0.483046, 0.493776, 0.504430, 0.515009, 0.525509, 0.535931, 0.546271, 0.556528, 0.566702, 0.576789, 0.586790, 0.596702, 0.606524, 0.616253, 0.625890, 0.635432, 0.644878, 0.654227, 0.663477, 0.672626, 0.681674, 0.690618, 0.699458, 0.708193, 0.716820, 0.725339, 0.733748, 0.742047, 0.750233, 0.758306, 0.766264, 0.774106, 0.781832, 0.789439, 0.796926, 0.804293, 0.811539, 0.818662, 0.825661, 0.832535, 0.839284, 0.845905, 0.852399, 0.858764, 0.864999, 0.871103, 0.877076, 0.882916, 0.888622, 0.894194, 0.899631, 0.904932, 0.910096, 0.915122, 0.920010, 0.924759, 0.929369, 0.933837, 0.938165, 0.942350, 0.946394, 0.950294, 0.954050, 0.957662, 0.961130, 0.964452, 0.967628, 0.970658, 0.973541, 0.976278, 0.978866, 0.981306, 0.983599, 0.985742, 0.987736, 0.989581, 0.991277, 0.992822, 0.994218, 0.995463, 0.996558, 0.997502, 0.998295, 0.998937, 0.999428, 0.999768, 0.999958, 0.999995, 0.999882, 0.999617, 0.999202, 0.998635, 0.997917, 0.997049, 0.996029, 0.994859, 0.993539, 0.992068, 0.990448, 0.988678, 0.986758, 0.984689, 0.982471, 0.980105, 0.977590, 0.974928, 0.972118, 0.969162, 0.966058, 0.962809, 0.959414, 0.955874, 0.952190, 0.948362, 0.944390, 0.940275, 0.936019, 0.931620, 0.927081, 0.922402, 0.917584, 0.912626, 0.907531, 0.902298, 0.896929, 0.891425, 0.885785, 0.880012, 0.874106, 0.868067, 0.861898, 0.855598, 0.849168, 0.842611, 0.835925, 0.829114, 0.822177, 0.815116, 0.807931, 0.800625, 0.793197, 0.785650, 0.777984, 0.770200, 0.762299, 0.754284, 0.746154, 0.737912, 0.729558, 0.721093, 0.712520, 0.703839, 0.695051, 0.686159, 0.677162, 0.668064, 0.658864, 0.649565, 0.640167, 0.630673, 0.621083, 0.611400, 0.601624, 0.591757, 0.581801, 0.571756, 0.561626, 0.551410, 0.541111, 0.530730, 0.520269, 0.509729, 0.499112, 0.488420, 0.477654, 0.466816, 0.455907, 0.444929, 0.433884, 0.422773, 0.411598, 0.400361, 0.389064, 0.377708, 0.366295, 0.354826, 0.343304, 0.331729, 0.320105, 0.308432, 0.296713, 0.284948, 0.273141, 0.261292, 0.249404, 0.237478, 0.225517, 0.213521, 0.201493, 0.189434, 0.177347, 0.165233, 0.153094, 0.140932, 0.128748, 0.116545, 0.104325, 0.092088, 0.079838, 0.067576, 0.055303, 0.043022, 0.030735, 0.018443, 0.006148, -0.006148, -0.018443, -0.030735, -0.043022, -0.055303, -0.067576, -0.079838, -0.092088, -0.104325, -0.116545, -0.128748, -0.140932, -0.153094, -0.165233, -0.177347, -0.189434, -0.201493, -0.213521, -0.225517, -0.237479, -0.249404, -0.261293, -0.273141, -0.284949, -0.296713, -0.308432, -0.320105, -0.331730, -0.343304, -0.354826, -0.366295, -0.377708, -0.389064, -0.400362, -0.411599, -0.422773, -0.433884, -0.444929, -0.455907, -0.466816, -0.477654, -0.488420, -0.499112, -0.509729, -0.520269, -0.530730, -0.541111, -0.551410, -0.561626, -0.571757, -0.581801, -0.591757, -0.601624, -0.611400, -0.621084, -0.630673, -0.640168, -0.649565, -0.658864, -0.668064, -0.677163, -0.686159, -0.695051, -0.703839, -0.712520, -0.721093, -0.729558, -0.737912, -0.746154, -0.754284, -0.762299, -0.770200, -0.777984, -0.785650, -0.793197, -0.800625, -0.807932, -0.815116, -0.822177, -0.829114, -0.835926, -0.842611, -0.849168, -0.855598, -0.861898, -0.868068, -0.874106, -0.880012, -0.885786, -0.891425, -0.896929, -0.902298, -0.907531, -0.912626, -0.917584, -0.922402, -0.927082, -0.931621, -0.936019, -0.940275, -0.944390, -0.948362, -0.952190, -0.955874, -0.959414, -0.962809, -0.966058, -0.969162, -0.972118, -0.974928, -0.977590, -0.980105, -0.982471, -0.984689, -0.986758, -0.988678, -0.990448, -0.992068, -0.993539, -0.994859, -0.996029, -0.997049, -0.997917, -0.998635, -0.999202, -0.999617, -0.999882, -0.999995, -0.999958, -0.999768, -0.999428, -0.998937, -0.998295, -0.997502, -0.996558, -0.995463, -0.994218, -0.992822, -0.991277, -0.989581, -0.987736, -0.985742, -0.983599, -0.981306, -0.978866, -0.976277, -0.973541, -0.970658, -0.967628, -0.964452, -0.961130, -0.957662, -0.954050, -0.950294, -0.946394, -0.942350, -0.938165, -0.933837, -0.929368, -0.924759, -0.920010, -0.915122, -0.910096, -0.904932, -0.899631, -0.894194, -0.888622, -0.882915, -0.877076, -0.871103, -0.864999, -0.858764, -0.852399, -0.845905, -0.839284, -0.832535, -0.825661, -0.818662, -0.811539, -0.804293, -0.796926, -0.789438, -0.781831, -0.774106, -0.766264, -0.758306, -0.750233, -0.742047, -0.733748, -0.725339, -0.716820, -0.708193, -0.699458, -0.690618, -0.681673, -0.672626, -0.663476, -0.654227, -0.644878, -0.635432, -0.625890, -0.616253, -0.606523, -0.596702, -0.586790, -0.576789, -0.566702, -0.556528, -0.546271, -0.535930, -0.525509, -0.515009, -0.504430, -0.493775, -0.483046, -0.472244, -0.461370, -0.450426, -0.439414, -0.428336, -0.417193, -0.405988, -0.394720, -0.383393, -0.372008, -0.360567, -0.349071, -0.337523, -0.325923, -0.314274, -0.302578, -0.290836, -0.279050, -0.267222, -0.255353, -0.243446, -0.231502, -0.219523, -0.207511, -0.195467, -0.183394, -0.171293, -0.159166, -0.147015, -0.134842, -0.122649, -0.110437, -0.098208, -0.085965, -0.073708, -0.061440, -0.049163, -0.036879, -0.024589, -0.012295, 0.000000};
