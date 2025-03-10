typedef struct mindlin {
double x, xprime,xdobleprime, k, b, c, f0, T, p0;
}Mindlin;

typedef struct gardner {
  double x, xprime,oldxprime,xdobleprime, K, Pb, a0, b0, t, M, K_scale, K_scalex, Pb_scalex, D, D2, Pb_scale, T, freq, ofreq;
}Gardner;

void init_mindlin(Mindlin* mind, double b, double k, double c);
void no_newsay();
int16_t mindlin_get_sample();

void init_gardner(Gardner* gd, int16_t kk, int16_t pb);
int16_t gardner_get_sample();

void init_balloon();
int16_t balloon_get_sample();
