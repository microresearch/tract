#include "sai.h"
#include "main.h"
#include <string.h>
#include <stdlib.h>
#include "wavetable.h"
#include "process.h"

//#define __UHSDR_DMAMEM __attribute__((section(".dma_buffer")))//ALIGN_32BYTES //__attribute__ ((section (".dmamem"))) 
#define ADC_CONVERTED_DATA_BUFFER_SIZE   ((uint32_t)  3)   /* Size of array aADCxConvertedData[] */

extern ADC_HandleTypeDef AdcHandle;
extern uint16_t aADCxConvertedData[ADC_CONVERTED_DATA_BUFFER_SIZE];
//extern uint16_t   aADCxConvertedData[ADC_CONVERTED_DATA_BUFFER_SIZE] __attribute__((section(".dma_buffer"))) __attribute__((aligned(0x20)));

/*
typedef struct {
    __packed audio_data_t l;
    __packed audio_data_t r;
} AudioSample_t;
*/

typedef struct
{
    AudioSample_t out[2*AUDIO_BLOCK_SIZE]; // 2 x 32 samples =64 left and 64 right
    AudioSample_t in[2*AUDIO_BLOCK_SIZE]; // now 48x2x2x4=48*16=768bytes
} dma_audio_buffer_t;

//static dma_audio_buffer_t dma __attribute__((section(".RAM_D3")));// __attribute__((aligned (32))); // was dma_mem
//static dma_audio_buffer_t dma __attribute__((section(".RAM_D2"))) __attribute__((aligned (4))); // NON//
//__attribute__((section(".dma_buffer"))) dma_audio_buffer_t dma;
static dma_audio_buffer_t dma __attribute__((section(".dma_buffer"))) __attribute__((aligned (32))); // was dma_mem
extern Wavetable wavtable;

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

__IO uint16_t vale = 0;

/* void UhsdrHwI2s_Codec_ClearTxDmaBuffer() */
/* { */
/*     memset((void*)&dma.iq_buf.out, 0, sizeof(dma.iq_buf.out)); */
/* } */

// #define PROFILE_APP
static void MchfHw_Codec_HandleBlock(uint16_t which)
{
	static float smoothedvalue;
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
	int16_t pbuf[2*AUDIO_BLOCK_SIZE], linbuf[2*AUDIO_BLOCK_SIZE],rinbuf[2*AUDIO_BLOCK_SIZE];
	static float test; int16_t testy;
    AudioSample_t *audio = &dma.in[offset];
    AudioSample_t *audioDst = &dma.out[offset];

    // Handle
	//    AudioDriver_I2SCallback(audio, iq, audioDst, sz);
	// void AudioDriver_I2SCallback(AudioSample_t *audio, IqSample_t *iq, AudioSample_t *audioDst, int16_t blockSize)
	// in audio_driver.c - test first with samples out or/???
	// question of stereo -> l and r (int16_t) so is audio->l, audio->r? for in and audioDst->l for out
	// 32 samples
	//int_to_floot

	// split left and right incoming

	// test ADC
  
	
	for (x=0;x<sz;x++){
		//		test+=2000;
		//		if (test>10000) test=-10000;
		//testy=audio[x].l;
		//		if (testy>10000) testy=10000;
		//		if (testy<-10000) testy=-10000;
		//		testy-=32000;
		linbuf[x]=testy;
		rinbuf[x]=audio[x].r;
	}
	//	int_to_floot(linbuf,flinbuffer,sz);
	//	int_to_floot(rinbuf,frinbuffer,sz);
	//	process(flinbuffer,frinbuffer,fbuffer,sz); // in process.c
	//dowavetable(fbuffer, &wavtable, 440.0f, sz); 
	//	floot_to_int(pbuf,fbuffer,sz);
	//	audio_comb_stereo(sz, audioDst, pbuf, pbuf);

	/// test aADCxConvertedData
	//	vol=(((float)(aADCxConvertedData[0]>>6)))/1024.0; // 16 bits -> 10 bits 
	//	vol=1.0;
	//	float value =(float)aADCxConvertedData[0]/65536.0f; 
	//	smoothedvalue += 0.1f * (value - smoothedvalue);
	//	int16_t vale=(int)(smoothedvalue*65536.0f);

ADC_ChannelConfTypeDef        sConfig;

 sConfig.Channel      = ADC_CHANNEL_9;//ADCx_CHANNEL;                /* Sampled channel number */
  sConfig.Rank         = ADC_REGULAR_RANK_1;          /* Rank of sampled channel number ADCx_CHANNEL */
  sConfig.SamplingTime = ADC_SAMPLETIME_8CYCLES_5;    /* Sampling time (number of clock cycles unit) */
  sConfig.SingleDiff   = ADC_SINGLE_ENDED;            /* Single-ended input channel */
  sConfig.OffsetNumber = ADC_OFFSET_NONE;             /* No offset subtraction */ 
  sConfig.Offset = 0;                                 /* Parameter discarded because offset correction is disabled */

  if (HAL_ADC_ConfigChannel(&AdcHandle, &sConfig) != HAL_OK)
  {
    /* Channel Configuration Error */
    Error_Handler();
  }

  
	/*##-3- Start the conversion process #######################################*/
  if (HAL_ADC_Start(&AdcHandle) != HAL_OK)
  {
    /* Start Conversation Error */
    Error_Handler();
  }

  if (HAL_ADC_PollForConversion(&AdcHandle, 10) != HAL_OK)
  {
    /* End Of Conversion flag not set on time */
    Error_Handler();
  }
  else
  {
  vale = HAL_ADC_GetValue(&AdcHandle);
  }
	vol=1.0;
	for (x=0;x<sz;x++){
		//		audioDst[x].l=(int)((float)linbuf[x]*vol);
		audioDst[x].l=vale>>1;//aADCxConvertedData[0];
		//		audioDst[x].r=((rand()%32768));
				audioDst[x].r=(int)((float)linbuf[x]*vol); // both l and r for the moment
	}
}

void HAL_SAI_RxCpltCallback(SAI_HandleTypeDef *hi2s)
{
	// these are the callbacks which call into handleblock and another handler...
	// but in the original they only handle A2 but this is maybe to sync it!
    MchfHw_Codec_HandleBlock(0);
}

void HAL_SAI_RxHalfCpltCallback(SAI_HandleTypeDef *hi2s)
{
    MchfHw_Codec_HandleBlock(1);
}

static void UhsdrHWI2s_Sai32Bits(SAI_HandleTypeDef* hsai)
{
    hsai->hdmarx->Init.PeriphDataAlignment = DMA_PDATAALIGN_WORD;
    hsai->hdmarx->Init.MemDataAlignment = DMA_MDATAALIGN_WORD;
    HAL_DMA_Init(hsai->hdmarx);

    HAL_SAI_InitProtocol(hsai, SAI_I2S_STANDARD, SAI_PROTOCOL_DATASIZE_32BIT, 2);
}

static void UhsdrHwI2s_SetBitWidth()
{
	UhsdrHWI2s_Sai32Bits(&hsai_BlockA1);
	UhsdrHWI2s_Sai32Bits(&hsai_BlockB1);
}



void UhsdrHwI2s_Codec_StartDMA()
{
	//		    UhsdrHwI2s_SetBitWidth();
	//SCB_CleanInvalidateDCache(); 
    // we clean the buffers since we don't know if we are in a "cleaned" memory segement
	memset((void*)&dma,0,sizeof(dma));

	//    memset((void*)&dma.iq_buf,0,sizeof(dma.iq_buf));
	HAL_SAI_Receive_DMA(&hsai_BlockA1,(int16_t*)dma.in,(sizeof(dma.in)/(sizeof(dma.in[0].l))));
	  //	  HAL_SAI_Receive_DMA(&hsai_BlockA1,(uint8_t*)dma.in,(sizeof(dma.in)/2));
	  //	  HAL_SAI_Receive_DMA(&hsai_BlockA1,dma.in,sizeof(dma.in)/sizeof(dma.in[0].l));
	//	SCB_CleanDCache_by_Addr((uint32_t*)(((uint32_t)dma.out) & ~(uint32_t)0x1F), sizeof(dma.out)/sizeof(dma.out[0].l +32));
	//SCB_InvalidateDCache_by_Addr((uint32_t *) &aADCxConvertedData[ADC_CONVERTED_DATA_BUFFER_SIZE/2], ADC_CONVERTED_DATA_BUFFER_SIZE);
	HAL_SAI_Transmit_DMA(&hsai_BlockB1,(int16_t*)dma.out,(sizeof(dma.out)/(sizeof(dma.out[0].l))));
		  //HAL_SAI_Transmit_DMA(&hsai_BlockB1,(uint8_t*)dma.out,(sizeof(dma.out)/2));
		  
		  //  HAL_SAI_Transmit_DMA(&hsai_BlockB1,dma.out,sizeof(dma.out)/sizeof(dma.out[0].l));
}


void UhsdrHwI2s_Codec_StopDMA(void)
{
    HAL_SAI_DMAStop(&hsai_BlockA1);
    HAL_SAI_DMAStop(&hsai_BlockB1);
}
