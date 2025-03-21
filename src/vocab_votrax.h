enum Phonemes {	EH3, EH2, EH1, PA0, DT, A1,  A2,  ZH, 	AH2,    I3,  I2,  I1,  M,  N,   B,   V,	CH,     SH,  Z,   AW1, NG, AH1, OO1, OO,	L,      K,   J,   H,   G,  F,   D,   S,	A,      AY,  Y1,  UH3, AH, P,   O,   I,	U,      Y,   T,   R,   E,  W,   AE,  AE1,	AW2,    UH2, UH1, UH,  O2, O1,  IU,  U1,	THV,    TH,  ER,  EH,  E1, AW,  PA1, STOP};

enum Intonations {I_1 =0x40,   I_2 =0x80,   I_3 =0xC0};

#define END 0xFF

// vocabularies - [0] is length: - get rid of that!

// Wizard of Wor

const unsigned char p1[] = {I_2|PA0, K, I, I3, L, PA0, W, O, O1, R, L, UH, K, PA0, F, O1, R, D, UH, B, UH3, L, S, K, O, O1, R, PA1, I_2|PA0, END}; //30,  

const unsigned char p2[] = {I, F, Y1, I3, IU, U, G, EH, EH3, T,
		      T, U, U1, P, AH1, W, ER, F, UH1, L, 
			PA1, I_2|PA0, AH1, EH3, I3, Y, L, T, A, K, 
			K, EH, R, UH1, V, Y1, IU, U, PA0, M, 
									  AH1, I3, Y, S, EH, L, F, PA1, I_2|PA0, END}; //50, 

const unsigned char p3[] = {I_2|PA0, Y1, IU, U, UH, R, I, N, PA1, I_2|PA0, END}; //11, 

const unsigned char p4[] = {I_2|PA0, THV, UH, PA0, D, UH, N, J, EH1, N, 
									  Z, UH, V, W, O, O1, R, R, PA1, I_2|PA1, END}; //21, 
const unsigned char p5[] = {AH1, UH3, I3, Y, AE1, EH3, M, PA1, END}; //9, 

const unsigned char p6[] = {PA1, THV, I_1|UH, W, I_1|I, Z, ER, D, I_1|UH, I_1|V, 
									  PA0, I_1|W, O, O1, R, R, PA1, END}; //18, 

const unsigned char p7[] = {W, UH, N, B, UH3, AH2, Y, T, F, R, UH, M, M, AH1, I3, Y1, P, R, I, T, Y, Y1, S, PA1, PA1, AE1, EH3, N, D, Y1, IU, U, L, EH, K, PA0, S, P, L, O, U1, D, PA1, END}; //44, 

const unsigned char p8[] = {AH1, EH3, I3, Y, K, R, E1, T, CH, ER, S, AH, UH3, R, R, A2, D, Y, O1, U1, AE1, EH3, K, T, I1, V, PA1, END}; //29,M, 

const unsigned char p9[] = {O, R, L, UH, K, PA0, W, I, L, PA0, EH, S, PA0, K, A2, I3, Y1, P, PA0, THV, R, U, THV, UH, D, O, O1, R, PA1 , END}; //32,W, 

const unsigned char p10[] = {IU, U, W, O1, N, T, H, I_1|AE1, I2, V, PA0, A2, T, CH, I_1|AE, N, S, PA0, F, O2, R, Y, O1, R, D, I_1|AE, N, S, PA1, END}; //31,Y1, 

const unsigned char p11[] = {R, E1, M, I_1|EH, M, B, ER, PA1, I_1|AH1, EH3, I3, Y, M, THV, UH, W, I, Z, ER, D, PA1, N, AH1, T, Y, IU, U1, U1, PA1, I_2|PA0, END}; //32,I_2|PA0, 

const unsigned char p12[] = {I_1|F, Y, I3, U1, K, AE1, EH3, N, T, PA1, B, I_1|E1, T, THV, UH1, R, EH, S, T, PA1, I_2|THV, I_1|EH1, N, Y, IU, U1, L, N, EH1, V, ER, PA1, G, I_1|EH1, T, THV, UH1, B, EH, S, T, PA0, I_2|PA0, END}; //45,I_1|I1, 

const unsigned char p13[] = {I, F, Y, I3, U, D, E1, S, T, R, O1, UH3, I3, AY, M, AH2, AH2, EH3, I3, Y, B, I_1|A, B, Y, Y1, S, PA1, AH1, I3, Y1, L, P, I_1|AH1, P, Y, IU, U, I1, N, THV, AY, AY, UH, V, EH, N, PA1, I_2|PA0, END}; //50,I_2|PA0, 

const unsigned char p14[] = {N, I_1|AH1, UH3, U1, AH1, EH3, I3, Y, M, G, EH, T, I, NG, M, AE, EH3, D, PA1, I_2|PA0, END}; //22,I_2|PA0, 

const unsigned char p15[] = {Y1, IU, U, W, I, L, N, I_1|EH, V, ER, L, E1, E1, V, W, I_1|O, I_1|O1, R, UH1, L, I_1|AH1, EH3, Y, V, PA1, I_2|PA0, END}; //28,I_2|PA0, 

const unsigned char p16[] = {G, I_1|AH1, R, W, O1, R, PA1, G, O1, O1, I_1|AE1, EH3, F, T, ER, THV, EH, M, PA1, I_2|PA0, END}; //22,I_2|PA0, 

const unsigned char p17[] = {W, AH1, T, CH, THV, UH, R, A, D, AH1, R, PA0, I_2|PA0, END}; //15,I_2|PA0, 

const unsigned char p18[] = {O, R, Y, ER, PA1, END}; //7,W, 

const unsigned char p19[] = {N, AH1, I_1|UH3, I_1|U1, Y1, IU, U1, U1, G, EH, T, PA0, TH, UH1, PA0, I_1|H, I_1|EH1, I_1|I3, V, E1, W, I_1|A1, I3, Y1, T, S, PA1, I_2|PA0, END}; //30,I_2|PA0, 

const unsigned char p20[] = {Y, O2, O2, R, I_1|AE1, EH3, I_1|S, I_1|K, I, NG, F, O, R, T, R, I_1|UH, I_1|B, UH3, L, PA1, I_2|PA0, END}; //23,I_2|PA0, 

const unsigned char p21[] = {F, Y, IU, U, T, R, AH1, EH3, I3, Y, I_1|AE1, N, Y, H, I_1|AH1, I_1|R, D, ER, PA1, I_2|PA0, Y1, IU, U, L, I_1|O1, N, L, Y, M, E, E1, T, W, I, TH, I_1|D, U, U, M, PA1, I_2|PA0, END}; //43,I, 

const unsigned char p22[] = {B, I_1|ER, I_1|R, W, O, R, PA1, G, I_1|AH1, I_1|R, W, O, R, PA1, AE1, EH3, N, D, TH, I_1|O, I_1|R, W, O, R, PA1, W, I, L, D, IU, I_1|U, Y1, IU, U, I, N, PA1, I_2|PA0, END}; //40,I_2|PA0, 

const unsigned char p23[] = {M, AH1, EH3, I3, Y, W, I_1|O, I_1|R, L, I, NG, S, AH1, UH3, R, V, I_1|EH, I_1|EH3, R, Y, V, I_1|EH, I_1|EH3, R, Y, H, I_1|UH, N, G, R, Y, PA1, I_2|PA0, END}; //35,I_2|PA0, 

const unsigned char p24[] = {AH1, I_1|EH3, I_1|I3, I_1|Y, M, I_1|AE1, D, J, I1, K, I, Z, S, T, R, I_1|AW, N, G, R, THV, AE1, EH3, N, Y, O2, O2, R, W, I_1|EH, I_1|P, UH1, N, S, PA1, END}; //36,M, 

const unsigned char p25[] = {Y, O2, O2, R, I_1|B, O, O2, I_1|N, S, W, I, L, PA0, L, AH1, EH3, Y, I, N, THV, UH, D, I_1|UH, N, J, EH1, N, S, UH, V, W, O, O1, R, R, PA1, I_2|PA0, END}; //39,I_2|PA0, 

const unsigned char p26[] = {AH1, EH3, I3, L, Y1, I_1|IU, I_1|U, D, E1, V, EH1, L, UH3, P, T, S, AH1, I3, AY, EH, N, S, PA1, W, E1, Y, D, E1, V, EH1, L, UH3, P, T, M, AE1, EH3, D, J, I1, K, PA1, END}; //44,W, 

const unsigned char p27[] = {I_1|A, I_1|I1, I_1|Y1, PA1, PA1, I, N, S, I_1|ER, I_1|T, PA1, I_1|K, I_1|O1, O2, I3, Y1, N, PA1, END}; //20,H, 

const unsigned char p28[] = {I_1|AH1, I_1|I3, I_1|Y, N, D, M, E, E1, PA1, END}; //11,F, 

const unsigned char p29[] = {I_1|I3, I_1|Y, M, PA0, AH2, O1, U1, D, AH1, PA0, S, P, AH2, I_1|I1, I_1|Y, T, PA1, END}; //19,AH1, 

const unsigned char p30[] = {EH, T, R, EH, D, Y, PA1, END}; //9,G, 

const unsigned char p31[] = {IU, U, D, PA0, B, I_1|EH1, T, ER, PA0, H, O, P, Y1, I_1|IU, I_1|U, D, O, N, T, I_1|F, I_1|AH1, I3, Y1, N, D, I_1|M, E, E1, PA1, PA1, PA1, END}; //33,Y1, 

const unsigned char p32[] = {N, UH, TH, ER, PA0, K, O1, O2, I3, Y1, N, F, O, R, M, I_1|AH1, I_1|I3, I_1|Y1, T, R, EH1, ZH, ER, T, CH, EH, S, T, PA1, END}; //31,AH1, 

const unsigned char p33[] = {H, I_1|AH1, H, I_1|AH1, H, AH1, H, AH1, PA1, END}; //11,PA1, 

const unsigned char p34[] = {AH2, PA0, I_1|G, I_1|IU, I_1|IU, IU, IU, D, PA1, M, AH1, I3, Y, P, I_1|EH, T, S, W, ER, R, G, EH, T, I, NG, H, I_1|UH, I_1|NG, G, R, Y, PA1, PA1, END}; //35,I_1|AH, 

const unsigned char p35[] = {Y1, IU, U, L, G, EH, T, PA1, PA1, THV, E, PA0, I_1|AH2, R, E, N, AH1, PA1, I_2|PA1, END}; //21,I_2|PA0, 

const unsigned char p36[] = {H, AH1, H, AH1, H, AH1, H, AH1, PA1, I_2|PA0, END}; //12,I_2|PA1, 

const unsigned char p37[] = {N, UH1, THV, ER, PA0, W, O, R, Y, ER, F, O1, R, M, AH1, I3, Y1, B, A, Y, B, Y1, Y, S, T, U, D, E, V, AH1, O2, U1, R, PA1, END}; //36,AH1, 

const unsigned char p38[] = {I_1|E, P, PA0, G, O, I1, Y1, NG, PA0, AE, N, D, Y, IU, U, W, I, L, F, I_1|AH1, I1, Y1, N, D, M, E, E1, PA1, END}; //30,K, 

const unsigned char p39[] = {F, E1, U, U, M, O, R, D, UH, N, J, EH1, N, S, PA1, AH1, N, D, Y, IU, I_1|U, I_1|L, B, E, E1, PA0, A, A2, END}; //30,AH1, 

const unsigned char p40[] = {I_1|O, I_1|R, L, O, R, D, PA1, END}; //9,W, 

const unsigned char p41[] = {W, I_1|O, I_1|R, L, O, R, D, PA1, I_2|PA0, END}; //11,I_2|PA0, 

const unsigned char p42[] = {AH1, M, B, AE, K, F, O, R, M, O, O1, R, PA1, PA1, W, I, TH, END}; //19,K, 

const unsigned char p43[] = {THV, UH, D, I_1|UH, I_1|N, J, EH, N, S, PA0, UH, V, PA0, W, O, O1, R, R, AH1, W, I_1|A2, I_1|AY, Y, T, PA0, Y, O, O1, R, R, I3, E1, T, I_1|ER, R, N, PA1, I_2|PA0, END}; //40,I_2|PA0, 

const unsigned char p44[] = {D, I_1|E, E1, P, I, N, THV, UH, PA0, K, AE, V, ER, N, S, PA0, UH, V, W, I_1|O, I_1|O1, R, R, PA1, PA1, Y1, IU, U, U1, W, I, L, M, I_1|E, I_1|E1, T, PA0, M, E, E1, PA0, I_2|PA0, END}; //44,I_2|PA0, 

const unsigned char p45[] = {TH, TH, AE1, EH3, NG, K, S, PA0, Y, IU, U, U1, PA1, END}; //15,PA1, 

const unsigned char p46[] = {Y1, IU, U, N, I_1|O1, I_1|O1, O1, Y1, IU, U, K, AE1, EH3, N, D, IU, U, B, I_1|EH, T, ER, PA1, I_2|PA0, END}; //25,I_2|PA0, 

const unsigned char p47[] = {I_1|ER, I_1|R, Y, B, AE1, EH3, K, PA1, PA1, AH1, EH3, I3, Y, K, AE1, EH3, N, T, W, I_1|A2, I_1|AY, Y, T, T, IU, U1, D, I_1|IU, U, I, T, UH1, G, I_1|A1, I_1|EH1, N, PA1, END}; //39,H, 

const unsigned char p48[] = {IU, U, K, AE1, EH3, N, S, T, I_1|AH1, R, T, AH1, N, I_1|IU, U1, W, PA1, I_2|PA1, B, UH, T, F, O1, R, N, AH1, I_1|UH3, I_1|U1, Y, O2, O2, R, TH, R, I_1|U1, U1, PA1, I_2|PA1, END}; //40,Y1, 

const unsigned char p49[] = {I_1|E, H, I_1|E, H, I_1|E, H, O, H, O, H, O, H, AH1, H, AH1, H, AH1, H, AH1, PA1, THV, I_1|AE, EH3, T, PA0, PA0, W, UH, Z, F, UH, N, PA1, END}; //35,H, 

const unsigned char p50[] = {I_1|EH, L, K, UH, M, PA1, T, U1, M, AH1, I1, Y, W, O1, I_1|ER, I_1|L, D, UH, V, W, O, O1, R, PA1, END}; //26,W, 

const unsigned char p51[] = {I_1|O, Y, IU, U1, V, K, I_1|UH, I_1|M, T, U1, S, K, O, O1, R, PA1, I1, N, THV, UH1, W, O1, I_1|ER, L, D, UH, V, W, I_1|O, O1, R, PA1, END}; //34,S, 

const unsigned char p52[] = {O2, O2, R, AW, F, T, IU, U1, I_2|S, E1, Y, THV, UH, I_2|W, I, Z, ER, D, PA1, THV, UH, M, AE1, EH3, D, J, I1, K, UH1, L, I_2|W, I, Z, ER, D, I_2|UH, V, W, O1, O2, R, PA1, PA1, END}; //45,Y, 

const unsigned char p53[] = {B, ER, R, W, O, R, H, AE, S, N, T, E, T, EH1, N, EH, N, Y, W, UH, N, PA0, I1, N, M, UH, N, TH, S, PA1, I_2|PA0, END}; //33,I_2|PA0, 

const unsigned char p54[] = {AH1, I3, Y, B, I_1|A, B, Y, Y1, S, PA0, B, R, E1, Y, TH, F, I_1|AH1, EH3, AY, R, PA1, END}; //23,M, 

const unsigned char p55[] = {AH1, EH3, I3, Y, L, F, R, AH1, I_1|EH3, I_1|EH3, I_1|Y, Y1, I3, U1, W, I1, TH, M, AH1, EH3, I3, Y, L, UH3, I_1|AH2, I_1|Y, T, N, I, NG, B, O, L, T, S, PA1, I_2|PA0, END}; //39,I_2|PA0, 

const unsigned char p56[] = {AH1, R, W, O, R, AE1, EH3, N, D, TH, O, R, W, O, R, PA1, B, Y, K, UH, M, PA1, I_2|PA1, I, N, V, I_1|I1, I_1|Z, I1, B, L, PA1, I_2|PA0, END}; //35,G, 

const unsigned char p57[] = {TH, O, R, W, O, R, PA1, I, Z, R, EH, D, PA1, M, E1, AY, N, PA1, AE1, EH3, N, D, H, I_1|UH, I_1|NG, G, R, Y, F, O, R, S, P, A2, I3, Y, S, F, U1, U1, D, PA1, I_2|PA0, END}; //45,I_2|PA0, 

const unsigned char p58[] = {O, R, Y, ER, F, AY, I2, R, PA1, AH1, EH3, I3, Y, D, R, AW, N, AY, I2, R, PA1, E1, T, CH, T, AH1, I3, Y1, M, AH1, EH3, I3, Y, UH1, P, I_1|AY, I_1|I3, R, PA1, END}; //41,W, 

const unsigned char p59[] = {O, R, Y, ER, PA1, I_2|PA0, END}; //8,I_2|W, 

const unsigned char p60[] = {IU, U1, V, PA0, J, UH, S, T, B, EH, N, F, R, AH1, I1, Y1, D, PA1, B, AH1, I2, Y1, END}; //24,Y, 

const unsigned char p61[] = {B, UH3, AH1, Y, T, THV, UH, B, O1, O1, L, T, PA1, I_2|PA0, END}; //16,I_2|PA0, 

const unsigned char p62[] = {UH, S, N, T, THV, AE1, T, L, UH3, I_1|AH2, I_1|Y, T, N, I, NG, B, O, L, T, PA0, D, E, L, I1, SH, UH1, S, PA1, END}; //30,W, 

const unsigned char p63[] = {AE, N, D, PA0, M, AH1, I1, Y1, T, EH1, L, EH1, P, O, R, T, I1, NG, S, P, EH, L, PA0, K, AE1, N, B, E, PA0, PA0, E1, V, EH, N, F, AE, S, T, ER, PA1, I_2|PA0, END}; //43,I_2|PA0, 

const unsigned char p64[] = {N, I_1|AH1, UH3, U1, Y, IU, U1, U1, N, O, O, THV, UH, T, A, A1, S, T, PA0, UH1, V, M, AH1, I2, Y1, M, AE1, EH3, D, J, I1, K, PA0, I_2|PA0, END}; //36,I_2|PA0, 

const unsigned char p65[] = {A, Y, B, E, Y, IU, U1, L, S, E, E1, M, E1, E, PA0, UH, G, A2, EH2, N, PA1, END}; //23,M, 

const unsigned char p66[] = {O2, O2, R, EH1, K, I_1|S, P, L, I_1|O, ZH, UH, N, PA0, W, UH1, S, M, I_1|Y1, I_1|U, S, I, K, T, U, M, AH1, I3, Y1, PA0, E1, ER, R, S, PA1, END}; //36,Y1, 

const unsigned char p67[] = {EH3, I3, Y, L, S, A, Y1, I, T, UH, G, A2, EH2, N, PA1, END}; //17,AH1, 

const unsigned char p68[] = {B, E1, E, F, O, R, W, O1, O, R, N, D, PA1, Y, IU, U, U1, PA1, UH1, P, R, O, O1, T, CH, PA1, THV, UH, P, I, T, PA1, I_2|PA0, END}; //35,I_2|PA0, 

const unsigned char p69[] = {Y, O2, O2, R, P, AE, TH, PA0, L, E, D, S, PA0, D, ER, EH1, K, T, L, Y1, PA0, T, U, U1, PA1, PA1, THV, UH, P, I, T, PA1, I_2|PA0, END}; //35,I_2|PA0, 

const unsigned char p70[] = {D, E1, E, P, ER, PA1, EH, V, ER, D, E1, E, P, ER, PA1, I, N, T, U, PA1, I_2|PA0, END}; //23,I_2|PA0, 

const unsigned char p71[] = {B, Y, W, EH, R, PA1, Y, IU, U, AH1, R, I, N, THV, UH, W, O, R, L, O, R, D, PA0, D, UH, N, J, EH1, N, Z, PA1, I_2|PA0, END}; //34,I_2|PA0, 

const unsigned char p72[] = {AH, AH1, PA1, Y, IU, U, TH, AW, T, Y, IU, U, K, OO, D, H, AH1, I2, Y1, D, PA1, B, UH, T, AH1, EH3, I3, Y, M, PA1, PA1, THV, UH, D, UH, N, J, EH1, N, M, AE, S, T, ER, PA1, I_2|PA0, END}; //48,I_2|PA0, 

const unsigned char p73[] = {TH, O, O1, R, PA0, B, ER, R, PA0, G, AH, R, PA1, D, I_1|I, I_1|N, ER, S, PA0, R, I_1|EH, I3, D, Y, PA1, I_2|PA0, END}; //28,I_2|PA0, 

const unsigned char p74[] = {I_1|A, I_1|I1, I_1|Y1, PA1, PA1, Y, O2, O2, R, S, P, A2, AY, Y, S, PA0, B, U, U1, T, S, PA0, UH, N, T, AH1, I2, Y1, D, PA1, END}; //32,H, 

const unsigned char p75[] = {M, AH1, EH3, I3, Y1, B, E, E1, S, T, S, PA0, R, UH, N, W, AH1, I2, Y1, L, D, PA1, I, N, THV, UH, W, O, R, L, O, R, D, PA1, D, UH, N, J, EH1, N, Z, PA1, I_2|PA0, END}; //45,I_2|PA0, 

const unsigned char p76[] = {AH1, I_1|UH3, I_1|U1, Y, O2, O2, R, O, N, L, Y, T, CH, AE, N, S, I1, S, Y, O2, O2, R, D, AE, N, S, PA1, END}; //29,N, 

const unsigned char p77[] = {AH, R, PA0, PA0, Y1, IU, U, PA0, I_2|PA0, F, I, T, PA0, PA0, T, U, PA0, I_2|PA0, S, ER, V, AH2, I2, Y1, V, PA0, PA0, THV, UH, PA0, I_2|PA0, P, I, T, PA1, END}; //37,I_2|PA0, 

const unsigned char p78[] = {P, S, PA1, PA1, AH1, UH3, I3, Y, M, UH, S, T, H, AE1, V, F, O, R, G, AH1, T, EH1, N, THV, UH, W, AW, L, S, PA1, END}; //32,U, 

const unsigned char p79[] = {W, AE1, ER, AH1, R, Y1, IU, U, G, O, I1, NG, T, U, H, AH1, I2, Y1, D, PA0, N, AH1, UH3, U, PA1, I_2|PA0, END}; //28,I_2|PA0, 
 
const unsigned char gorf0[] = {6, 33, 34, END };
const unsigned char gorf1[] = {6, 60, END };
const unsigned char gorf2[] = {50, 28, 2, 9, 13, END };
const unsigned char gorf3[] = {47, 1, 12, END };
const unsigned char gorf4[] = {46, 12, END };
const unsigned char gorf5[] = {47, 0, 13, 30, END };
const unsigned char gorf6[] = {50, 13, 8, 11, 41, 24, 6, 17, 49, 13, END };
const unsigned char gorf7[] = {8, 13, 51, 57, 58, END };
const unsigned char gorf8[] = {21, 13, 51, 56, 43, 43, END };
const unsigned char gorf9[] = {21, 43, END };
const unsigned char gorf10[] = {51, 43, END };
const unsigned char gorf11[] = {50, 15, 2, 0, 13, 13, 30, 26, 58, END };
const unsigned char gorf12[] = {14, 46, 0, 30, END };
const unsigned char gorf13[] = {14, 46, 2, 30, END };
const unsigned char gorf14[] = {14, 44, END };
const unsigned char gorf15[] = {14, 0, 2, 13, END };
const unsigned char gorf16[] = {14, 50, 4, 42, 1, 13, 13, END };
const unsigned char gorf17[] = {25, 51, 30, 1, 42, END };
const unsigned char gorf18[] = {25, 47, 13, 50, 42, END };
const unsigned char gorf19[] = {25, 47, 0, 37, 42, 9, 13, END };
const unsigned char gorf20[] = {25, 43, 8, 35, 13, 9, 25, 35, 52, 24, 18, END };
const unsigned char gorf21[] = {25, 53, 35, 9, 60, 13, END };
const unsigned char gorf22[] = {25, 53, 35, 9, 60, 13, 31, END };
const unsigned char gorf23[] = {25, 58, 13, 8, 24, END };
const unsigned char gorf24[] = {25, 21, 13, 16, 59, 31, 13, 59, 31, END };
const unsigned char gorf25[] = {30, 60, 29, 2, 13, 30, 58, END };
const unsigned char gorf26[] = {30, 60, 31, 42, 43, 53, 11, 41, END };
const unsigned char gorf27[] = {30, 60, 31, 42, 43, 53, 11, 34, 30, END };
const unsigned char gorf28[] = {30, 55, 40, 12, END };
const unsigned char gorf29[] = {30, 43, 19, 31, END };
const unsigned char gorf30[] = {59, 12, 37, 8, 9, 41, 43, END };
const unsigned char gorf31[] = {59, 13, 30, END };
const unsigned char gorf32[] = {2, 13, 2, 12, 41, END };
const unsigned char gorf33[] = {2, 31, 25, 5, 60, 37, END };
const unsigned char gorf34[] = {29, 24, 59, 28, 17, 39, 37, END };
const unsigned char gorf35[] = {29, 38, 43, END };
const unsigned char gorf36[] = {28, 50, 24, 46, 25, 42, 11, 25, END };
const unsigned char gorf37[] = {30, 26, 1, 13, 58, 35, 24, END };
const unsigned char gorf38[] = {28, 4, 53, 53, 43, 43, 29, 29, END };
const unsigned char gorf39[] = {28, 4, 53, 43, 29, 41, 6, 13, END };
const unsigned char gorf40[] = {28, 4, 53, 43, 29, 60, 1, 13, END };
const unsigned char gorf41[] = {28, 4, 53, 43, 29, 41, 6, 13, 31, END };
const unsigned char gorf42[] = {27, 21, 27, 21, 27, 21, 27, 21, END };
const unsigned char gorf43[] = {27, 51, 43, 30, 58, END };
const unsigned char gorf44[] = {27, 47, 0, 15, END };
const unsigned char gorf45[] = {27, 11, 42, 11, 20, END };
const unsigned char gorf46[] = {21, 11, 41, END };
const unsigned char gorf47[] = {21, 11, 34, END };
const unsigned char gorf48[] = {11, 12, 37, 21, 31, 11, 14, 24, END };
const unsigned char gorf49[] = {39, 13, 0, END };
const unsigned char gorf50[] = {39, 13, 31, 58, 42, END };
const unsigned char gorf51[] = {11, 31, END };
const unsigned char gorf52[] = {24, 11, 15, END };
const unsigned char gorf53[] = {24, 61, 20, END };
const unsigned char gorf54[] = {12, 44, 60, 42, END };
const unsigned char gorf55[] = {12, 40, 55, 15, END };
const unsigned char gorf56[] = {12, 8, 11, 41, END };
const unsigned char gorf57[] = {12, 8, 9, 41, END };
const unsigned char gorf58[] = {13, 39, 43, END };
const unsigned char gorf59[] = {13, 59, 25, 31, 42, END };
const unsigned char gorf60[] = {13, 35, 8, 41, 31, END };
const unsigned char gorf61[] = {13, 38, END };
const unsigned char gorf62[] = {13, 21, 53, 55, 45, END };
const unsigned char gorf63[] = {62, END };
const unsigned char gorf64[] = {37, 24, 5, 33, 58, END };
const unsigned char gorf65[] = {37, 43, 60, 37, 47, 58, END };
const unsigned char gorf66[] = {37, 43, 11, 31, 9, 13, 0, 43, 31, END };
const unsigned char gorf67[] = {37, 43, 49, 12, 38, 42, 0, 30, END };
const unsigned char gorf68[] = {37, 22, 54, 17, END };
const unsigned char gorf69[] = {43, 53, 55, 14, 21, 42, END };
const unsigned char gorf70[] = {43, 53, 55, 14, 21, 42, 31, END };
const unsigned char gorf71[] = {43, 53, 55, 14, 21, 35, 42, 31, END };
const unsigned char gorf72[] = {31, 44, 25, END };
const unsigned char gorf73[] = {17, 39, 37, END };
const unsigned char gorf74[] = {17, 21, 35, 42, END };
const unsigned char gorf75[] = {31, 51, 12, END };
const unsigned char gorf76[] = {31, 37, 6, 9, 41, 31, END };
const unsigned char gorf77[] = {3, END };
const unsigned char gorf78[] = {31, 58, 15, 8, 11, 34, 15, 35, 24, END };
const unsigned char gorf79[] = {42, 32, 25, END };
const unsigned char gorf80[] = {56, 51, END };
const unsigned char gorf81[] = {56, 50, END };
const unsigned char gorf82[] = {57, 51, END };
const unsigned char gorf83[] = {42, 21, 0, 41, 12, END };
const unsigned char gorf84[] = {42, 40, END };
const unsigned char gorf85[] = {42, 54, 55, END };
const unsigned char gorf86[] = {49, 13, 14, 60, 33, 42, 35, 14, 35, 24, END };
const unsigned char gorf87[] = {45, 38, 43, 33, 34, 0, 43, END };
const unsigned char gorf88[] = {45, 38, 43, 33, 34, 0, 43, 31, END };
const unsigned char gorf89[] = {45, 11, 24, END };
const unsigned char gorf90[] = {34, 9, 55, END };
const unsigned char gorf91[] = {41, 54, 55, 55, END };
const unsigned char gorf92[] = {41, 11, 55, 55, END };
const unsigned char gorf93[] = {34, 54, 55, 55, END };
const unsigned char gorf94[] = {34, 11, 55, 55, END };
const unsigned char gorf95[] = {41, 38, 43, END };
const unsigned char gorf96[] = {41, 55, 53, 43, 31, 2, 24, 29, END };
const unsigned char gorf97[] = {29, 53, 43, END };
const unsigned char gorf98[] = {29, 52, 43, END };
const unsigned char gorf99[] = {45, 39, 24, END };
const unsigned char gorf100[] = {28, 4, 53, 43, 15, 41, 6, 13, END };
const unsigned char gorf101[] = {25, 53, 35, 9, 33, 13, 13, END };
const unsigned char gorf102[] = {50, 42, 46, 0, 25, END };
const unsigned char gorf103[] = {14, 8, 9, 34, 42, END };
const unsigned char gorf104[] = {25, 21, 13, 25, 58, END };
const unsigned char gorf105[] = {30, 41, 15, 21, 55, 58, END };
const unsigned char gorf106[] = {30, 51, 31, 42, END };
const unsigned char gorf107[] = {28, 47, 24, 50, 25, 31, 41, END };
const unsigned char gorf108[] = {28, 21, 0, 42, END };
const unsigned char gorf109[] = {37, 21, 11, 43, END };
const unsigned char gorf110[] = {42, 43, 8, 11, 41, END };
const unsigned char gorf111[] = {31, 55, 37, 43, 60, 13, END };
const unsigned char gorf112[] = {61, 24, END };
const unsigned char gorf113[] = {27, 6, 33, 24, END };
const unsigned char gorf114[] = {2, 12, 37, 21, 11, 43, END };

// list of vocabs

const unsigned char *vocablist_wow[79] =  {p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79};

const unsigned char *vocablist_gorf[115] =  {gorf0,gorf1,gorf2,gorf3,gorf4,gorf5,gorf6,gorf7,gorf8,gorf9,gorf10,gorf11,gorf12,gorf13,gorf14,gorf15,gorf16,gorf17,gorf18,gorf19,gorf20,gorf21,gorf22,gorf23,gorf24,gorf25,gorf26,gorf27,gorf28,gorf29,gorf30,gorf31,gorf32,gorf33,gorf34,gorf35,gorf36,gorf37,gorf38,gorf39,gorf40,gorf41,gorf42,gorf43,gorf44,gorf45,gorf46,gorf47,gorf48,gorf49,gorf50,gorf51,gorf52,gorf53,gorf54,gorf55,gorf56,gorf57,gorf58,gorf59,gorf60,gorf61,gorf62,gorf63,gorf64,gorf65,gorf66,gorf67,gorf68,gorf69,gorf70,gorf71,gorf72,gorf73,gorf74,gorf75,gorf76,gorf77,gorf78,gorf79,gorf80,gorf81,gorf82,gorf83,gorf84,gorf85,gorf86,gorf87,gorf88,gorf89,gorf90,gorf91,gorf92,gorf93,gorf94,gorf95,gorf96,gorf97,gorf98,gorf99,gorf100,gorf101,gorf102,gorf103,gorf104,gorf105,gorf106,gorf107,gorf108,gorf109,gorf110,gorf111,gorf112,gorf113,gorf114};
