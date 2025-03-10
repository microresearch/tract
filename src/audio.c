#include "stm32h7xx_hal.h"
#include "main.h"
#include <string.h>
#include <stdlib.h>
#include "wavetable.h"
#include "samplerate.h"
#include "resources.h"
#include "tms5200x.h"
#include "sp0256.h"
#include "vot.h"
#include "process.h"
#include "klatt_phoneme.h"
#include "parwave.h"
#include "holmes.h"
#include "nvp.h"
#include "rs.h"
#include "wormed.h"
#include "tube.h"
#include "raven.h"
#include "LF.h"

extern DMA_HandleTypeDef hdma_sai1_a;
extern SAI_HandleTypeDef hsai_BlockA1;
extern ADC_HandleTypeDef AdcHandle;

#define DMA_MAX_SZE                     0xFFFF
#define DMA_MAX(_X_)                (((_X_) <= DMA_MAX_SZE)? (_X_):DMA_MAX_SZE)
#define AUDIODATA_SIZE                  2   /* 16-bits audio data size */

AudioSample_t audio_out[2*AUDIO_BLOCK_SIZE] __attribute__((section(".dma_buffer"))) __attribute__((aligned (32))); // was dma_mem; // 2 x 32 samples =64 left and 64 right or 2x48
extern Wavetable wavtable;

__IO uint16_t adc_buffer[3];

int16_t signall[128];

float smoothed_adc_value[5]={0.0f, 1.0f, 0.0f, 0.0f, 0.0f}; // as mode is inverted and we don't want to start with COMPOSTF?

float _mode, _speed, _selx=0.5, _sely=1.0, _selz=0.5;

typedef struct wormer_ {
  u8 maxextent;
  float sampleratio;
  int16_t(*getsample)(void);
  void(*newsay)(void);
  u8 xy;
  u8 TTS;
} wormer;

static const wormer tmser={0, 1.0f, tms_get_sample, tms_newsay, 0, 0};
static const wormer tmsphoner={0, 1.0f, tms_get_sample_allphon, tms_newsay_allphon, 0, 0};
static const wormer sp0256er={0, 0.325f, sp0256_get_sample, sp0256_newsay, 0, 0};
static const wormer sp0256er1={0, 0.325f, sp0256_get_samplevocabbankone, sp0256_newsayvocabbankone, 0, 0};
//sp0256_get_samplevocabbankone

static const wormer votraxer={0, 1.0f, votrax_get_sample, votrax_newsay, 0, 0};
static const wormer votraxGer={0, 1.0f, votrax_get_samplegorf, votrax_newsaygorf, 0, 0};


//rsynth_2005
static const wormer rsynthy={0, 1.0f, rsynth_get_sample, rsynth_newsay, 0, 0};
static const wormer rsynthy1={0, 1.0f, rsynth_get_sample_single, rsynth_newsay_single, 0, 0};

//rsynth-2.0-port
static const wormer klatter={0, 1.0f, klatt_get_sample, klatt_newsay, 0, 0};  // elements
//static const wormer simpleklatter={38, 0.5f, simpleklatt_get_sample, simpleklatt_newsay, 2, 0}; 
static const wormer klattsingle={0, 1.0f, klatt_get_sample_single, klatt_newsay_single, 0, 0};
static const wormer klattvocab={0, 1.0f, klatt_get_sample_vocab, klatt_newsay_vocab, 0, 0};

//nvp.c
static const wormer nvper={0, 0.5f, nvp_get_sample, nvp_newsay, 0, 0};
static const wormer nvperv={0, 0.5f, nvp_get_sample_vocab, nvp_newsay_vocab, 0, 0};

static const wormer simpleklatter={38, 0.5f, simpleklatt_get_sample, simpleklatt_newsay, 2, 0}; 

// do wavetable
static const wormer waver={0, 0.2f, wave_get_sample1, wave_newsay, 0, 0};
static const wormer waver1={0, 0.2f, wave_get_sample2, wave_newsay, 0, 0};
static const wormer waver2={0, 0.2f, wave_get_sample3, wave_newsay, 0, 0};
static const wormer waver3={0, 0.2f, wave_get_sample4, wave_newsay, 0, 0};
static const wormer waver4={0, 0.2f, wave_get_sample5, wave_newsay, 0, 0};
static const wormer waver5={0, 0.2f, wave_get_sample6, wave_newsay, 0, 0};

// wormed.c sources
static const wormer chirper={0, 0.4f, chirp_get_sample, wormed_newsay, 0, 0};
static const wormer votgloter={0, 0.4f, votglot_get_sample, wormed_newsay, 0, 0};
static const wormer impulser={0, 0.8f, impulsive_get_sample, wormed_newsay, 0, 0};
static const wormer smper={0, 0.1f, sampled_source_getsample, wormed_newsay, 0, 0}; // not using 2nd adc _selz
static const wormer nater={0, 0.8f, natural_source_getsample, wormed_newsay, 0, 0};
static const wormer impulsz={0, 0.8f, impulsivex_get_sample, wormed_newsay, 0, 0};

// tube.c
static const wormer tuber={0, 0.8f, tube_get_sample, tube_newsay, 0, 0};

//raven.c mindlin
static const wormer mindlin={0, 1.0f, mindlin_get_sample, no_newsay, 0, 0};
static const wormer gardner={0, 2.0f, gardner_get_sample, no_newsay, 0, 0};
static const wormer balloner={0, 1.0f, balloon_get_sample, no_newsay, 0, 0};

// lfgen and lfgen2.c
static const wormer LFF={0, 1.0f, LF_get_sample, no_newsay, 0, 0};
static const wormer rosen={0, 1.0f, rosenberg_get_sample, no_newsay, 0, 0};
static const wormer shimmer={0, 0.4f, shimmer_get_sample, no_newsay, 0, 0};
static const wormer klglottt={0, 1.0f, klg_get_sample, no_newsay, 0, 0};
static const wormer LF2={0, 1.0f, LF2_get_sample, no_newsay, 0, 0};

static const wormer *wormlist[]={&tmsphoner, &tmser, &votraxer, &votraxGer, &sp0256er, &sp0256er1, &rsynthy, &rsynthy1, &klatter, &nvper, &nvperv, &waver, &waver1, &waver2, &waver3, &waver4, &waver5, &chirper, &votgloter, &impulser, &smper, &nater, &impulsz, &tuber, &mindlin, &gardner, &balloner, &LFF, &rosen, &shimmer, &klglottt, &LF2};  // 32

/*
- DONEadd LF(test first), 
- DONEadditional wavetables: 6,7,8,9,10
- DONEothers from stripped_worm:  
- test each and check params
- add in mode adc
- program and test all
*/

#define MAXED(var, max) \
  if (var > (max)) {	\
    var = (max);	\
  }

#define CONSTRAIN(var, min, max) \
  if (var < (min)) { \
    var = (min); \
  } else if (var > (max)) { \
    var = (max); \
  }

static inline void doadc(){
  float value;

  ADC_ChannelConfTypeDef        sConfig;
	sConfig.Channel      = ADC_CHANNEL_16;//ADCx_CHANNEL;                /* Sampled channel number */
	sConfig.Rank         = ADC_REGULAR_RANK_1;          /* Rank of sampled channel number ADCx_CHANNEL */
	sConfig.SamplingTime = ADC_SAMPLETIME_8CYCLES_5;    /* Sampling time (number of clock cycles unit) */
	sConfig.SingleDiff   = ADC_SINGLE_ENDED;            /* Single-ended input channel */
	sConfig.OffsetNumber = ADC_OFFSET_NONE;             /* No offset subtraction */ 
	sConfig.Offset = 0;                                 /* Parameter discarded because offset correction is disabled */

  if (HAL_ADC_ConfigChannel(&AdcHandle, &sConfig) != HAL_OK)
  {
    /* Channel Configuration Error */
    //    Error_Handler();
  }

  
  if (HAL_ADC_Start(&AdcHandle) != HAL_OK)
  {
    /* Start Conversation Error */
    //    Error_Handler();
  }

  if (HAL_ADC_PollForConversion(&AdcHandle, 10) != HAL_OK)
  {
    /* End Of Conversion flag not set on time */
    //    Error_Handler();
  }
  else
  {
  adc_buffer[0] = HAL_ADC_GetValue(&AdcHandle);
  }

  	sConfig.Channel      = ADC_CHANNEL_17;//ADCx_CHANNEL;                /* Sampled channel number */
	sConfig.Rank         = ADC_REGULAR_RANK_1;          /* Rank of sampled channel number ADCx_CHANNEL */
	sConfig.SamplingTime = ADC_SAMPLETIME_8CYCLES_5;    /* Sampling time (number of clock cycles unit) */
	sConfig.SingleDiff   = ADC_SINGLE_ENDED;            /* Single-ended input channel */
	sConfig.OffsetNumber = ADC_OFFSET_NONE;             /* No offset subtraction */ 
	sConfig.Offset = 0;                                 /* Parameter discarded because offset correction is disabled */

  if (HAL_ADC_ConfigChannel(&AdcHandle, &sConfig) != HAL_OK)
  {
    /* Channel Configuration Error */
    //    Error_Handler();
  }

  
  if (HAL_ADC_Start(&AdcHandle) != HAL_OK)
  {
    /* Start Conversation Error */
    //    Error_Handler();
  }

  if (HAL_ADC_PollForConversion(&AdcHandle, 10) != HAL_OK)
  {
    /* End Of Conversion flag not set on time */
    //    Error_Handler();
  }
  else
  {
  adc_buffer[1] = HAL_ADC_GetValue(&AdcHandle);
  }
	sConfig.Channel      = ADC_CHANNEL_14;//ADCx_CHANNEL;                /* Sampled channel number */
	sConfig.Rank         = ADC_REGULAR_RANK_1;          /* Rank of sampled channel number ADCx_CHANNEL */
	sConfig.SamplingTime = ADC_SAMPLETIME_8CYCLES_5;    /* Sampling time (number of clock cycles unit) */
	sConfig.SingleDiff   = ADC_SINGLE_ENDED;            /* Single-ended input channel */
	sConfig.OffsetNumber = ADC_OFFSET_NONE;             /* No offset subtraction */ 
	sConfig.Offset = 0;                                 /* Parameter discarded because offset correction is disabled */

  if (HAL_ADC_ConfigChannel(&AdcHandle, &sConfig) != HAL_OK)
  {
    /* Channel Configuration Error */
    //    Error_Handler();
  }

  
  if (HAL_ADC_Start(&AdcHandle) != HAL_OK)
  {
    /* Start Conversation Error */
    //    Error_Handler();
  }

  if (HAL_ADC_PollForConversion(&AdcHandle, 10) != HAL_OK)
  {
    /* End Of Conversion flag not set on time */
    //    Error_Handler();
  }
  else
  {
  adc_buffer[2] = HAL_ADC_GetValue(&AdcHandle);
  }
  
  value =(float)adc_buffer[0]/65536.0f; 
  smoothed_adc_value[0] += 0.1f * (value - smoothed_adc_value[0]);
  _speed=smoothed_adc_value[0];
  CONSTRAIN(_speed,0.0f,1.0f);

  value =(float)adc_buffer[1]/65536.0f; 
  smoothed_adc_value[1] += 0.1f * (value - smoothed_adc_value[1]);
  _selz=smoothed_adc_value[1];
  CONSTRAIN(_selz,0.0f,1.0f);

  value =(float)adc_buffer[2]/65536.0f; 
  smoothed_adc_value[2] += 0.01f * (value - smoothed_adc_value[2]); // smoother
  _mode=smoothed_adc_value[2];
  CONSTRAIN(_mode,0.0f,1.0f);
}

static int16_t	monobuffer[256];


static inline void floot_to_int(int16_t* outbuffer, float* inbuffer,u16 howmany){
	int16_t tmp;

  for (int n = 0; n < howmany; n++) {
    tmp = inbuffer[n] * 32768.0f;
    tmp = (tmp <= -32768) ? -32768 : (tmp >= 32767) ? 32767 : tmp;
    outbuffer[n] = (int16_t)tmp;
		}
}

static inline void int_to_floot(int16_t* inbuffer, float* outbuffer, u16 howmany){
  for (int n = 0; n < howmany; n++) {
    outbuffer[n]=(float)(inbuffer[n])/32768.0f;
  }
}

static void Codec_HandleBlock(uint16_t which)
{
	int16_t sample;
	float samplespeed;
	float value;
	u16 samplespeedref;

	// with callback to audio processing handle to check out
	///	SCB_CleanDCache_by_Addr((uint32_t*)(((uint32_t)dma.out) & ~(uint32_t)0x1F), sizeof(dma.out)/sizeof(dma.out[0].l +32));
	//	SCB_InvalidateDCache_by_Addr((uint32_t*)(((uint32_t)dma.in) & ~(uint32_t)0x1F), sizeof(dma.in)/sizeof(dma.in[0].l +32));

	// Transfer complete interrupt
    // Point to 2nd half of buffers
	const size_t sz = IQ_BLOCK_SIZE;
	const uint16_t offset = which == 0?sz:0;
	uint16_t x, val;
	float fbuffer[2*AUDIO_BLOCK_SIZE],flinbuffer[2*AUDIO_BLOCK_SIZE],frinbuffer[2*AUDIO_BLOCK_SIZE];
	float vol;
	static u8 _intmode=0;
	u8 oldmode=1;
	int16_t pbuf[2*AUDIO_BLOCK_SIZE], linbuf[2*AUDIO_BLOCK_SIZE],rinbuf[2*AUDIO_BLOCK_SIZE];
	static float test; int16_t testy;
	AudioSample_t *audioDst = &audio_out[offset];

	doadc();
	samplespeedref=_speed*1028.0f;
	MAXED(samplespeedref, 1023);
	samplespeed=logspeed[samplespeedref];  
	
	//	vol=(float)smoothed_adc_value[0]/65536.0;
	// make into a call
	//	for (x=0;x<sz;x++){
	uint16_t triggered=0; // 16*,17,18,19,20,21DONE to test 6/3/2025 // TESTY
	  //	  	  audioDst[x].l=(int)((float)((rand()%32768))*vol);
	oldmode=_intmode;
	_intmode=(int)(_mode*35.0);
	if (_intmode>31) _intmode=31;
	if (oldmode!=_intmode) {triggered=1;};
	//	_intmode=31; triggered=0; //test each in turn TESTY
	samplerate_simple(monobuffer, samplespeed, sz, wormlist[_intmode]->getsample, wormlist[_intmode]->newsay , wormlist[_intmode]->sampleratio, triggered);
	  //		audioDst[x].l=vale>>1;//aADCxConvertedData[0];
	  //	  audioDst[x].l=((rand()%32768)); // left is our DAC out
	  //	  audioDst[x].r=((rand()%32768));
	  //			audioDst[x].r=(int)((float)linbuf[x]*vol); // both l and r for the moment
	  //	}
	  for (x=0;x<sz;x++){
	    audioDst[x].l=monobuffer[x];
	  }
}


void HAL_SAI_TxHalfCpltCallback(SAI_HandleTypeDef *hsai_BlockA1)
{
	if(hsai_BlockA1->Instance==SAI1_Block_A)
	{
//		HAL_GPIO_TogglePin(LD3_GPIO_Port, LD3_Pin);
//		HalfTransfer_CallBack_FS();
	  Codec_HandleBlock(1);
	}
}

//void HAL_I2S_TxCpltCallback(I2S_HandleTypeDef *hi2s)
//{
//	TransferComplete_CallBack_FS();
//}

void HAL_SAI_TxCpltCallback(SAI_HandleTypeDef *hsai_BlockA1)
{
	if(hsai_BlockA1->Instance==SAI1_Block_A)
	{
    Codec_HandleBlock(0);
//		HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
//		TransferComplete_CallBack_FS();
	}
}

void Codec_StartDMA()
{
	//		    UhsdrHwI2s_SetBitWidth();
	//SCB_CleanInvalidateDCache(); 
    // we clean the buffers since we don't know if we are in a "cleaned" memory segement
  memset((void*)&audio_out,0,sizeof(audio_out)); // ???
  HAL_SAI_Transmit_DMA(&hsai_BlockA1,(int16_t*)audio_out,(sizeof(audio_out)/(sizeof(audio_out[0].l))));
    
}
