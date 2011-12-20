//OBSOLUTE DEFINE
//#define X_SIZE	256
//#define Y_SIZE	256
//OBSOLUTE DEFINE END

//DEFINE VERSION
#define CORE_NAME	"SORA_CORE"
#define CORE_DEV_STATUS	"A.C.S PRE"
//STATUS IS NOT RELATED TO SVN/GOOGLE CODE
#define CORE_VERSION	"0.1.9a"
//VERSION IS NOT RELATED TO SVN/GOOGLE CODE
#define CORE_CODE_NAME	"KIRINO"
//CODE NAME IS NOT RELATED TO SVN/GOOGLE CODE
#define CORE_PLATFORM	"LINUX/UNIX"
//DEFINE VERSION END

//DEFINE DEBUG
//#define debug(msg)	printf("%s\n", msg) //ON
#define debug(msg)	{} //OFF
//DEFINE DEBUG END

//STANDARD MACRO
//#define tostring(var)do {         printf("%s: %s\n", #var, var);     } while (0)
//STANDARD MACRO END

//DEFINE CONSTANT
#define TIME  (time((time_t *) NULL))
#define X_EXP	8
#define Y_EXP	8
#define HIGH	255
#define LOW		0
#define OFFSET 128
#define TEXT_BUF_LEN 10000
#define DIV 8
#define XS (X_SIZE / DIV)
#define YS (Y_SIZE / DIV)
#define DTH 0.7
#define L_BASE 100
#define ROOT2 (double)1.41421356
//DEFINE CONSTANT END

