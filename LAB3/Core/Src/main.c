/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2022 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#define _USE_MATH_DEFINES
#define ARM_MATH_CM4
#include "arm_math.h"					//include library for matrix calculation
#include "math.h"						//include library for Pi value
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
TIM_HandleTypeDef htim1;
TIM_HandleTypeDef htim2;
TIM_HandleTypeDef htim3;

UART_HandleTypeDef huart2;

/* USER CODE BEGIN PV */
uint64_t _micros = 0;							//time in microsecond
float DegRel =0;								//Degree Relative
float DegAbs =0;								//Degree Absolute
float RadRel =0;								//Radian Relative
float RadAbs =0;								//Radian Absolute
float unwrap=0;									//Degree Unwrap
float unwraprad = 0;							//Radian Unwrap
//Variable for Unwrapping function
float Po ;
float Polast;
float Pnlast;
//Variable for velocity approximate function
float omega=0;
float tetalast=0;
//Variable for kalman filter
float Rad_predict = 0;
float Omega_predict = 0;
//sample time in kalman
const float delta_T = 0.01;

//Kalman All Matrix define

float32_t X_f32[3] =//initial condition
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  0.0,
  0.0,
  0.0,
};
float32_t U_f32[1] =//initial condition
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  0.0
};
float32_t A_f32[9] =
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  1.0,     delta_T,   delta_T*delta_T/2.0,
  0.0,     1.0,      delta_T,
  0.0,     0.0,      1.0

};
float32_t B_f32[3] =
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  0,
  0,
  0
};
float32_t G_f32[3] =
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  delta_T*delta_T*delta_T/6.0,
  delta_T*delta_T/2.0,
  delta_T
};
float32_t Q_f32[1] =
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  0.00828
};
//mesurement
float32_t C_f32[3] =
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  1,	0,		0
};
/*const float32_t D_f32[1] =
{
  // Const,   numTaps,   blockSize,   numTaps*blockSize
  0
};*/
float32_t R_f32[1] =
{
  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
  0.01
};
float32_t Z_f32[1] =//input unwrap
{
  //Const,   numTaps,   blockSize,   numTaps*blockSize
  0.0
};
float32_t Y_f32[1] =
{
  //Const,   numTaps,   blockSize,   numTaps*blockSize
  0.0
};
float32_t AT_f32[9];//3X3
float32_t CT_f32[3];//3X1
float32_t GT_f32[3];//1X3
float32_t Xp_f32[3];//3X1
float32_t Pp_f32[9];//3X3
float32_t P_f32[9]=
	{
		  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
		  0.0,	0.0,	0.0,
		  0.0,	0.0,	0.0,
		  0.0,	0.0,	0.0,
	};//3X3
float32_t I_f32[9]=
	{
		  /* Const,   numTaps,   blockSize,   numTaps*blockSize */
		  1.0,	0.0,	0.0,
		  0.0,	1.0,	0.0,
		  0.0,	0.0,	1.0,
	};//3X3
float32_t S_f32[1];//1X1
float32_t K_f32[3];//3X1
float32_t GQGt_f32[9];//3X3
float32_t GQ_f32[3];//3X1
//calculation
float32_t AX_f32[3];//3X1
float32_t BU_f32[3];//3X1
//Xp
float32_t AP_f32[9];//3X3
float32_t APAt_f32[9];//3X3
//Pp
float32_t AX_f32[3];//3X1
float32_t CXp_f32[1];//1X1
float32_t CPp_f32[3];//1X3
float32_t CPpCt_f32[3];//1X3

float32_t PpCt_f32[3];//3X1
float32_t S_inv_f32[1];//1X1
float32_t KY_f32[3];//3X1
float32_t KC_f32[9];//3X3
float32_t I_KC_f32[9];//3X3

arm_matrix_instance_f32 A;      /* Matrix A Instance */
arm_matrix_instance_f32 AT;     /* Matrix AT(A transpose) */
arm_matrix_instance_f32 B;      /* Matrix B Instance */
arm_matrix_instance_f32 G;      /* Matrix G Instance */
arm_matrix_instance_f32 GT;     /* Matrix GT(G transpose) */
arm_matrix_instance_f32 Q;      /* Matrix Q Instance */
arm_matrix_instance_f32 C;      /* Matrix C Instance */
arm_matrix_instance_f32 CT;     /* Matrix CT(C transpose) */
arm_matrix_instance_f32 R;      /* Matrix R Instance */
arm_matrix_instance_f32 Y;      /* Matrix Y Instance */
arm_matrix_instance_f32 X;      /* Matrix X Instance */
arm_matrix_instance_f32 Xp;     /* Matrix Xp(X predict) */
arm_matrix_instance_f32 P;      /* Matrix P Instance */
arm_matrix_instance_f32 Pp;     /* Matrix Pp(P predict) */
arm_matrix_instance_f32 U;      /* Matrix U Instance */
arm_matrix_instance_f32 K;      /* Matrix K(kalman gain) Instance */
arm_matrix_instance_f32 S;      /* Matrix S Instance */
arm_matrix_instance_f32 Z;      /* Matrix Z Instance */
arm_matrix_instance_f32 GQ;     /* Matrix GQ Instance */
arm_matrix_instance_f32 GQGt;   /* Matrix GQGt Instance */
arm_matrix_instance_f32 AX;     /* Matrix Ax Instance */
arm_matrix_instance_f32 BU;     /* Matrix BU Instance */
arm_matrix_instance_f32 APAt;   /* Matrix APAt Instance */
arm_matrix_instance_f32 AP;     /* Matrix AP Instance */
arm_matrix_instance_f32 CXp;    /* Matrix CXp Instance */
arm_matrix_instance_f32 CPp;    /* Matrix CPp Instance */
arm_matrix_instance_f32 CPpCt;  /* Matrix CPpCt Instance */
arm_matrix_instance_f32 PpCt;   /* Matrix PpCt Instance */
arm_matrix_instance_f32 S_inv;  /* Matrix S_inv Instance */
arm_matrix_instance_f32 KY;     /* Matrix KY Instance */
arm_matrix_instance_f32 I;      /* Matrix I Instance */
arm_matrix_instance_f32 KC;     /* Matrix I Instance */
arm_matrix_instance_f32 I_KC;   /* Matrix I Instance */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART2_UART_Init(void);
static void MX_TIM1_Init(void);
static void MX_TIM2_Init(void);
static void MX_TIM3_Init(void);
/* USER CODE BEGIN PFP */
//function prototype
uint64_t micros();
void UpdatePosition();
float unwapping(float Pn);
void velocityapproximate(float );
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_USART2_UART_Init();
  MX_TIM1_Init();
  MX_TIM2_Init();
  MX_TIM3_Init();
  /* USER CODE BEGIN 2 */

  uint32_t srcRows, srcColumns;
  //set rows and columns for matrix
  srcRows = 3;
  srcColumns = 3;
  arm_mat_init_f32(&A, srcRows, srcColumns, A_f32);
  arm_mat_init_f32(&AT, srcRows, srcColumns, AT_f32);
  arm_mat_trans_f32(&A, &AT);
  arm_mat_init_f32(&GQGt, srcRows, srcColumns, GQGt_f32);
  arm_mat_init_f32(&P, srcRows, srcColumns, P_f32);
  arm_mat_init_f32(&Pp, srcRows, srcColumns, Pp_f32);
  arm_mat_init_f32(&AP, srcRows, srcColumns, AP_f32);
  arm_mat_init_f32(&APAt, srcRows, srcColumns, APAt_f32);
  arm_mat_init_f32(&I, srcRows, srcColumns, I_f32);
  arm_mat_init_f32(&KC, srcRows, srcColumns, KC_f32);
  arm_mat_init_f32(&I_KC, srcRows, srcColumns, I_KC_f32);
  srcRows = 1;
  srcColumns = 3;
  arm_mat_init_f32(&C, srcRows, srcColumns, C_f32);
  arm_mat_init_f32(&GT, srcRows, srcColumns, GT_f32);
  arm_mat_init_f32(&CPp, srcRows, srcColumns, CPp_f32);
  arm_mat_init_f32(&CPpCt, srcRows, srcColumns, CPpCt_f32);
  srcRows = 3;
  srcColumns = 1;
  arm_mat_init_f32(&X, srcRows, srcColumns, X_f32);
  arm_mat_init_f32(&Xp, srcRows, srcColumns, Xp_f32);
  arm_mat_init_f32(&CT, srcRows, srcColumns, CT_f32);
  arm_mat_trans_f32(&C, &CT);
  arm_mat_init_f32(&G, srcRows, srcColumns, G_f32);
  arm_mat_trans_f32(&G, &GT);
  arm_mat_init_f32(&B, srcRows, srcColumns, B_f32);
  arm_mat_init_f32(&K, srcRows, srcColumns, K_f32);
  arm_mat_init_f32(&GQ, srcRows, srcColumns, GQ_f32);
  arm_mat_init_f32(&AX, srcRows, srcColumns, AX_f32);
  arm_mat_init_f32(&BU, srcRows, srcColumns, BU_f32);
  arm_mat_init_f32(&CPpCt, srcRows, srcColumns, CPpCt_f32);
  arm_mat_init_f32(&PpCt, srcRows, srcColumns, PpCt_f32);
  arm_mat_init_f32(&KY, srcRows, srcColumns, KY_f32);

  srcRows = 1;
  srcColumns = 1;
  arm_mat_init_f32(&Q, srcRows, srcColumns, Q_f32);
  arm_mat_init_f32(&U, srcRows, srcColumns, U_f32);
  arm_mat_init_f32(&R, srcRows, srcColumns, R_f32);
  arm_mat_init_f32(&Y, srcRows, srcColumns, Y_f32);
  arm_mat_init_f32(&S, srcRows, srcColumns, S_f32);
  arm_mat_init_f32(&Z, srcRows, srcColumns, Z_f32);
  arm_mat_init_f32(&S_inv, srcRows, srcColumns, S_inv_f32);
  arm_mat_mult_f32(&G,&Q,&GQ);
  arm_mat_mult_f32(&GQ,&GT,&GQGt);
  arm_mat_init_f32(&CXp, srcRows, srcColumns, CXp_f32);

  	HAL_TIM_Base_Start_IT(&htim3);
	HAL_TIM_Base_Start_IT(&htim2);
	HAL_TIM_Encoder_Start(&htim1, TIM_CHANNEL_1);
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
	  UpdatePosition();  //Update Degree and radian from QEI
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);
  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 8;
  RCC_OscInitStruct.PLL.PLLN = 100;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 4;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_3) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief TIM1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM1_Init(void)
{

  /* USER CODE BEGIN TIM1_Init 0 */

  /* USER CODE END TIM1_Init 0 */

  TIM_Encoder_InitTypeDef sConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM1_Init 1 */

  /* USER CODE END TIM1_Init 1 */
  htim1.Instance = TIM1;
  htim1.Init.Prescaler = 0;
  htim1.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim1.Init.Period = 64511;
  htim1.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim1.Init.RepetitionCounter = 0;
  htim1.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  sConfig.EncoderMode = TIM_ENCODERMODE_TI12;
  sConfig.IC1Polarity = TIM_ICPOLARITY_RISING;
  sConfig.IC1Selection = TIM_ICSELECTION_DIRECTTI;
  sConfig.IC1Prescaler = TIM_ICPSC_DIV1;
  sConfig.IC1Filter = 0;
  sConfig.IC2Polarity = TIM_ICPOLARITY_RISING;
  sConfig.IC2Selection = TIM_ICSELECTION_DIRECTTI;
  sConfig.IC2Prescaler = TIM_ICPSC_DIV1;
  sConfig.IC2Filter = 0;
  if (HAL_TIM_Encoder_Init(&htim1, &sConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim1, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM1_Init 2 */

  /* USER CODE END TIM1_Init 2 */

}

/**
  * @brief TIM2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM2_Init(void)
{

  /* USER CODE BEGIN TIM2_Init 0 */

  /* USER CODE END TIM2_Init 0 */

  TIM_ClockConfigTypeDef sClockSourceConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM2_Init 1 */

  /* USER CODE END TIM2_Init 1 */
  htim2.Instance = TIM2;
  htim2.Init.Prescaler = 99;
  htim2.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim2.Init.Period = 4294967295;
  htim2.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim2.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim2) != HAL_OK)
  {
    Error_Handler();
  }
  sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
  if (HAL_TIM_ConfigClockSource(&htim2, &sClockSourceConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim2, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM2_Init 2 */

  /* USER CODE END TIM2_Init 2 */

}

/**
  * @brief TIM3 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM3_Init(void)
{

  /* USER CODE BEGIN TIM3_Init 0 */

  /* USER CODE END TIM3_Init 0 */

  TIM_ClockConfigTypeDef sClockSourceConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM3_Init 1 */

  /* USER CODE END TIM3_Init 1 */
  htim3.Instance = TIM3;
  htim3.Init.Prescaler = 99;
  htim3.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim3.Init.Period = 9999;
  htim3.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim3.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim3) != HAL_OK)
  {
    Error_Handler();
  }
  sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
  if (HAL_TIM_ConfigClockSource(&htim3, &sClockSourceConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim3, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM3_Init 2 */

  /* USER CODE END TIM3_Init 2 */

}

/**
  * @brief USART2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART2_UART_Init(void)
{

  /* USER CODE BEGIN USART2_Init 0 */

  /* USER CODE END USART2_Init 0 */

  /* USER CODE BEGIN USART2_Init 1 */

  /* USER CODE END USART2_Init 1 */
  huart2.Instance = USART2;
  huart2.Init.BaudRate = 115200;
  huart2.Init.WordLength = UART_WORDLENGTH_8B;
  huart2.Init.StopBits = UART_STOPBITS_1;
  huart2.Init.Parity = UART_PARITY_NONE;
  huart2.Init.Mode = UART_MODE_TX_RX;
  huart2.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart2.Init.OverSampling = UART_OVERSAMPLING_16;
  if (HAL_UART_Init(&huart2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART2_Init 2 */

  /* USER CODE END USART2_Init 2 */

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : B1_Pin */
  GPIO_InitStruct.Pin = B1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(B1_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LD2_Pin */
  GPIO_InitStruct.Pin = LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LD2_GPIO_Port, &GPIO_InitStruct);

}

/* USER CODE BEGIN 4 */
void UpdatePosition(){
	static uint64_t TimestampCalculationLoop  =0;
	if(micros() -TimestampCalculationLoop >=1000){	//Update every 1000 Hz
		uint16_t RawRead = TIM1->CNT;
		DegRel = (RawRead/3072.0)*360;
		DegAbs =((RawRead%3072)/3072.0)*360;

		RadRel =(RawRead/3072.0)*(2*M_PI);
		RadAbs =((RawRead%3072)/3072.0)*(2*M_PI);

	}
}
void velocityapproximate(float teta){
	omega = 100.0*(teta-tetalast);
	tetalast = teta;
}
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim)
{
	if (htim == &htim2)
	{
		_micros += 4294967295;
	}
	if (htim == &htim3){
			unwrap = unwapping(DegRel);
			unwraprad = unwrap*M_PI/180.0;
			velocityapproximate(unwraprad);
			// Kalman Calcultation
			Z_f32[0] = unwraprad;
			//Xp
			arm_mat_mult_f32(&A,&X,&AX);
			arm_mat_mult_f32(&B,&U,&BU);
			arm_mat_add_f32(&AX,&BU,&Xp);
			//Pp
			arm_mat_mult_f32(&A,&P,&AP);
			arm_mat_mult_f32(&AP,&AT,&APAt);
			arm_mat_add_f32(&APAt,&GQGt,&Pp);
			//Y
			arm_mat_mult_f32(&C,&Xp,&CXp);
			arm_mat_sub_f32(&Z,&CXp,&Y);
			//S
			arm_mat_mult_f32(&C,&Pp,&CPp);
			arm_mat_mult_f32(&CPp,&CT,&CPpCt);
			arm_mat_add_f32(&CPpCt,&R,&S);
			//K
			arm_mat_mult_f32(&Pp,&CT,&PpCt);
			arm_mat_inverse_f32(&S,&S_inv);
			arm_mat_mult_f32(&PpCt,&S_inv,&K);
			//X
			arm_mat_mult_f32(&K,&Y,&KY);
			arm_mat_add_f32(&KY,&Xp,&X);
			//P
			arm_mat_mult_f32(&K,&C,&KC);
			arm_mat_sub_f32(&I,&KC,&I_KC);
			arm_mat_mult_f32(&I_KC,&Pp,&P);

			Rad_predict = X_f32[0];
			Omega_predict = X_f32[1];
		}
}
uint64_t micros()
{
	return _micros + htim2.Instance->CNT;
}
float unwapping(float Pn){		//Upwrapping
	static float Pmax = 7560;   //360*21 (21 rounds 1 unwrap)
	static float threshold = 7560*0.5;
	if (Pn-Pnlast <= -threshold){
		Po = Polast + Pmax;
	}
	else if (Pn - Pnlast >= threshold){
		Po = Polast - Pmax;
	}
	else{
		Po = Polast;
	}
	Polast = Po;
	Pnlast = Pn;
	return Pn+Po;
}
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */

