#define ADC_CONVERTED_DATA_BUFFER_SIZE   ((uint32_t)  32)   /* Size of array aADCxConvertedData[] */

//gcc test.c -otest -lm -std=gnu99 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#define randi() rand()

void main(void)
{
  int xx;
  xx=ADC_CONVERTED_DATA_BUFFER_SIZE;
  printf("xxx: %d\n",xx);
}
