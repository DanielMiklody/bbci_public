/*
  bbci_acquire_tmsi.c
  
  1. state = bbci_acquire_tmsi('init', state); 
  2. [data] = bbci_acquire_tmsi(state);
  3. [data, marker_time] = bbci_acquire_tmsi(state);
  4. [data, marker_time, marker_descr,state] = bbci_acquire_tmsi(state);
  5. bbci_acquire_tmsi('close');
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include <winsock2.h>
#include <winbase.h>
#include <windows.h>
#include <wchar.h>
#include <conio.h>
#include <tchar.h> 
#include <assert.h>
#include <queue>
#include <time.h>
#include <string>
#include <sstream>
using namespace std;
using std::string;
#pragma comment(lib,"ws2_32.lib")


#include "mex.h"
#include "TmsiSDK.h" 

#define NUM_CHAN 34
/*
	GLOBAL VARIABLES
*/
static bool g_bIsConnected = false;  // shows where we are connected to the server or not
static bool g_bIsMac = false; // shows whether its a direct connection
static SOCKET g_Socket; // active socket, which is connected to the server
static SOCKET g_ServerSocket; // active socket, which is connected to the server
static SOCKET g_ActiveSocket; 
//static FILE *g_Fp;

static int g_NumberOfChannels;
unsigned int g_BytesPerSample;
static ULONG g_SampleRateInMilliHz;
static ULONG g_SignalBufferSizeInSamples;
static ULONG g_SampleRateInHz;
static HANDLE g_Handle;
static HINSTANCE g_LibHandle;
static unsigned long long g_total;

// Afterwards the user will be able to set them(we will have to change the acquire function a little bit
const int p    = 5;          // Duration of the pause before acquisition (seconds, min 1s)
const int N    = 8;          // number of channels
const int bps  = 4;          // bytes per sample
const int Fs   = 500;        // sampling rate
const int Ns   = 20;         // number of samples to read in the buffer each time (between 20 and 250)
const int buff = bps*Ns*N;   // Number of bytes to read in total for all channels
int  datenPuffer[Ns*N]={0};

// the status name constants
static const char* FIELD_IP = "hostIP";
static const char* FIELD_MAC = "hostMAC";
static const char* FIELD_Channels = "numChan";
static const char* FIELD_Freq = "fs";


struct chData
{
  int channels[34];
  unsigned long long timeStamp;
};
struct markerData
{
	string description;
	unsigned long long pos;
	unsigned long long timeStamp;
	double time;
};

void quitTMSI();
static queue<chData> gDataQueue;
static queue<markerData> gMarkerQueue;

static HANDLE ghMutexData; 
static HANDLE ghMutexMarkers; 
static HANDLE ghMutexCount; 
static HANDLE  hServerThread;
static HANDLE  hObtainThread;
static unsigned long long g_CurCount=0;
static unsigned int g_numCh=8;
static bool bTerminate=false;

static POPEN fpOpen;
static PCLOSE fpClose; 
static PSTART fpStart;
static PSTOP fpStop;	
static PSETSIGNALBUFFER fpSetSignalBuffer;
static PGETSAMPLES	fpGetSamples;
static PGETSIGNALFORMAT fpGetSignalFormat; 
static PFREE fpFree;
static PLIBRARYINIT fpLibraryInit;
static PLIBRARYEXIT fpLibraryExit;
static PGETFRONTENDINFO fpGetFrontEndInfo;
static PSETRTCTIME fpSetRtcTime;
static PGETRTCTIME fpGetRtcTime;
static PSETRTCALARMTIME fpSetRtcAlarmTime;
static PGETRTCALARMTIME fpGetRtcAlarmTime;
static PGETERRORCODE fpGetErrorCode;
static PGETERRORCODEMESSAGE fpGetErrorCodeMessage;
static PFREEDEVICELIST fpFreeDeviceList;
static PGETDEVICELIST fpGetDeviceList;
static PGETCONNECTIONPROPERTIES fpGetConnectionProperties;
static PSETMEASURINGMODE fpSetMeasuringMode;
static PSETREFCALCULATION fpSetRefCalculation;
static PGETBUFFERINFO fpGetBufferInfo;

// Functions for Mobita
static PSTARTCARDFILE fpStartCardFile;
static PSTOPCARDFILE fpStopCardFile;
static PGETCARDFILESAMPLES fpGetCardFileSamples;
static PGETCARDFILESIGNALFORMAT fpGetCardFileSignalFormat;
static POPENCARDFILE fpOpenCardFile;
static PGETCARDFILELIST fpGetCardFileList;
static PCLOSECARDFILE fpCloseCardFile;
static PGETRECORDINGCONFIGURATION fpGetRecordingConfiguration;
static PSETRECORDINGCONFIGURATION fpSetRecordingConfiguration;
static PGETEXTFRONTENDINFO fpGetExtFrontEndInfo;

//Functions for Nexus10-MKII
static PGETRANDOMKEY fpGetRandomKey;
static PUNLOCKFRONTEND fpUnlockFrontEnd;
static PGETOEMSIZE fpGetOEMSize;
static PSETOEMDATA fpSetOEMData;
static PGETOEMDATA fpGetOEMData;
static PSETSTORAGEMODE fpSetStorageMode;

static POPENFIRSTDEVICE fpOpenFirstDevice;


static PSIGNAL_FORMAT psf = NULL;
static FRONTENDINFO FrontEndInfo;

static time_t g_startTime;

time_t getUnixTimeStamp()
{
	time_t seconds;
	FILETIME* ft = new FILETIME;
	SYSTEMTIME* st = new SYSTEMTIME;
	seconds = time (NULL);
	seconds*=1000;
	GetSystemTimeAsFileTime(ft);
	FileTimeToSystemTime(ft,st);
	seconds+=st->wMilliseconds;
	return seconds;
}

DWORD WINAPI threadObtain( LPVOID lpParam ) 
{

	unsigned int* SignalBuffer, SignalBufferSizeInBytes;
	float Fval[NUM_CHAN];
	SignalBufferSizeInBytes = g_SignalBufferSizeInSamples*NUM_CHAN*sizeof(SignalBuffer[0]);
	SignalBuffer = (unsigned int*) malloc(SignalBufferSizeInBytes);

			while(1) {
//			mexPrintf("bla");
			int BytesReturned = fpGetSamples(g_Handle,(PULONG) SignalBuffer, SignalBufferSizeInBytes);
			int NrSamples = BytesReturned/(g_NumberOfChannels*sizeof(unsigned int));
			for(int i=0;i<NrSamples;i++) {
				for(int j=0;j<g_NumberOfChannels;j++) {
					int ind = i*g_NumberOfChannels+j;
					if(SignalBuffer[ind] == OVERFLOW_32BITS && 
							(psf[j].Type == CHANNELTYPE_EXG || 
							psf[j].Type == CHANNELTYPE_BIP || 
							psf[j].Type == CHANNELTYPE_AUX ))
								Fval[j] = 0;
					else {
						switch(psf[j].Format) {
							case SF_UNSIGNED:
								Fval[j] = SignalBuffer[ind] *  psf[j].UnitGain +  psf[j].UnitOffSet ;
								break ;
							case SF_INTEGER: // signed integer
								Fval[j] = ((int) SignalBuffer[ind]) *  psf[j].UnitGain +  psf[j].UnitOffSet ;
								break ;
							default : 
								Fval[j] = 0 ; // For unknown types, set the value to zero 
								break ;
						}
					}
				}
				
				
				if(bTerminate)
				{
					bTerminate=false;
					return 0;
				}
				WaitForSingleObject(ghMutexData, INFINITE );
				WaitForSingleObject(ghMutexCount, INFINITE );
				
				chData newData;
				  for(int i=0;i<g_numCh;i++)
						newData.channels[i] = Fval[i];
				  newData.timeStamp = g_CurCount++;

					  gDataQueue.push(newData);
	
				
				ReleaseMutex(ghMutexData);
				ReleaseMutex(ghMutexCount);
			
				if(bTerminate)
				{
					bTerminate=false;
					return 0;
				}
			}
		}
	
	return 0;
}

char **DeviceList = NULL;
int NrOfDevices=0;

// initialize tmsi
void initTMSI() {
	TCHAR LibraryName[255] = _T("\\TmsiSDK.dll");
	TCHAR Path[	MAX_PATH ];
	int ErrorCode=0;
	GetSystemDirectory(Path, sizeof(Path) / sizeof(TCHAR) );
	lstrcat(Path, LibraryName);
	g_LibHandle = LoadLibrary(Path);
	

	if(!g_LibHandle) {
		mexPrintf("ERROR. Cannot load DLL. Are the tmsi drivers installed?");
		return;
	}
	
	fpOpen				= (POPEN)			GetProcAddress(g_LibHandle,"Open");
	fpClose				= (PCLOSE)			GetProcAddress(g_LibHandle,"Close");
	fpStart				= (PSTART)			GetProcAddress(g_LibHandle,"Start");
	fpStop				= (PSTOP)			GetProcAddress(g_LibHandle,"Stop");
	fpSetSignalBuffer	= (PSETSIGNALBUFFER)GetProcAddress(g_LibHandle,"SetSignalBuffer");
	fpGetSamples		= (PGETSAMPLES)		GetProcAddress(g_LibHandle,"GetSamples");
	fpGetBufferInfo		= (PGETBUFFERINFO)	GetProcAddress(g_LibHandle,"GetBufferInfo");
	fpGetSignalFormat	= (PGETSIGNALFORMAT)GetProcAddress(g_LibHandle,"GetSignalFormat"); 
	fpFree				= (PFREE)			GetProcAddress(g_LibHandle, "Free" ); 
	fpLibraryInit		= (PLIBRARYINIT)	GetProcAddress(g_LibHandle, "LibraryInit" ); 
	fpLibraryExit		= (PLIBRARYEXIT)	GetProcAddress(g_LibHandle, "LibraryExit" ); 
	fpGetFrontEndInfo	= (PGETFRONTENDINFO) GetProcAddress(g_LibHandle, "GetFrontEndInfo" ); 
	fpSetRtcTime		= (PSETRTCTIME)		GetProcAddress(g_LibHandle, "SetRtcTime" ); 
	fpGetRtcTime		= (PGETRTCTIME)		GetProcAddress(g_LibHandle, "GetRtcTime" ); 
	fpSetRtcAlarmTime	= (PSETRTCALARMTIME)GetProcAddress(g_LibHandle, "SetRtcAlarmTime" ); 
	fpGetRtcAlarmTime	= (PGETRTCALARMTIME)GetProcAddress(g_LibHandle, "GetRtcAlarmTime" ); 
	fpGetErrorCode		= (PGETERRORCODE)	GetProcAddress(g_LibHandle, "GetErrorCode" ); 
	fpGetErrorCodeMessage = (PGETERRORCODEMESSAGE) GetProcAddress(g_LibHandle, "GetErrorCodeMessage" ); 
	fpGetDeviceList		= (PGETDEVICELIST)	GetProcAddress(g_LibHandle, "GetDeviceList" ); 
	fpFreeDeviceList	= (PFREEDEVICELIST)	GetProcAddress(g_LibHandle, "FreeDeviceList" ); 
	fpStartCardFile		= (PSTARTCARDFILE)	GetProcAddress(g_LibHandle, "StartCardFile" ); 
	fpStopCardFile		= (PSTOPCARDFILE)	GetProcAddress(g_LibHandle, "StopCardFile" ); 
	fpGetCardFileSamples	= (PGETCARDFILESAMPLES)	GetProcAddress(g_LibHandle, "GetCardFileSamples" ); 
	fpGetConnectionProperties = (PGETCONNECTIONPROPERTIES)	GetProcAddress(g_LibHandle, "GetConnectionProperties" ); 
	fpGetCardFileSignalFormat = (PGETCARDFILESIGNALFORMAT) GetProcAddress(g_LibHandle, "GetCardFileSignalFormat" ); 
	fpOpenCardFile		= (POPENCARDFILE) GetProcAddress(g_LibHandle, "OpenCardFile" ); 
	fpGetCardFileList	= (PGETCARDFILELIST) GetProcAddress(g_LibHandle, "GetCardFileList" ); 
	fpCloseCardFile		= (PCLOSECARDFILE) GetProcAddress(g_LibHandle, "CloseCardFile" );
	fpGetExtFrontEndInfo = (PGETEXTFRONTENDINFO) GetProcAddress(g_LibHandle, "GetExtFrontEndInfo");
	fpSetMeasuringMode	= (PSETMEASURINGMODE) GetProcAddress(g_LibHandle, "SetMeasuringMode" );
	fpGetRecordingConfiguration = (PGETRECORDINGCONFIGURATION) GetProcAddress(g_LibHandle, "GetRecordingConfiguration" );
	fpSetRecordingConfiguration = (PSETRECORDINGCONFIGURATION) GetProcAddress(g_LibHandle, "SetRecordingConfiguration" );
	fpSetRefCalculation = (PSETREFCALCULATION) GetProcAddress(g_LibHandle, "SetRefCalculation" );
	fpGetRandomKey = (PGETRANDOMKEY) GetProcAddress(g_LibHandle, "GetRandomKey");
	fpUnlockFrontEnd=(PUNLOCKFRONTEND) GetProcAddress(g_LibHandle, "UnlockFrontEnd");
	fpGetOEMSize=(PGETOEMSIZE) GetProcAddress(g_LibHandle, "GetOEMSize");
	fpGetOEMData=(PGETOEMDATA) GetProcAddress(g_LibHandle, "GetOEMData");
	fpSetOEMData=(PSETOEMDATA) GetProcAddress(g_LibHandle, "SetOEMData");
	fpOpenFirstDevice = (POPENFIRSTDEVICE) GetProcAddress(g_LibHandle, "OpenFirstDevice" );
	fpSetStorageMode = (PSETSTORAGEMODE) GetProcAddress(g_LibHandle, "SetStorageMode");

	if(!fpGetRecordingConfiguration) {
		mexPrintf("Failed to load functions from dll");
		return;
	}
	
	g_Handle = fpLibraryInit( TMSiConnectionUSB, &ErrorCode );
	
	if(ErrorCode) {
		mexPrintf("Failed to initialize the library with Library Init. Errorcode = %d", ErrorCode); 
		return; 
	}

	DeviceList = fpGetDeviceList( g_Handle, &NrOfDevices);
	
		if(!NrOfDevices) {
			mexPrintf("0 devices found. Have you connected any devices?");
			fpLibraryExit( g_Handle );
			return;
		}

		char FrontEndName[MAX_FRONTENDNAME_LENGTH];
		
		psf = NULL;
		FRONTENDINFO FrontEndInfo;
		int Status;
		char *DeviceLocator = DeviceList[0];
		Status = fpOpen(g_Handle,DeviceLocator);
		if(!Status) {
			mexPrintf("Could not Open");
		}
		Status = fpGetFrontEndInfo(g_Handle,&FrontEndInfo);
		psf = fpGetSignalFormat(g_Handle,FrontEndName);
//		for(int i=0;i<34;i++)
//			mexPrintf("Offset ch %d, %f, %f",i,psf[i].UnitGain,psf[i].UnitOffSet);
			
		if(!psf) {
			ErrorCode = fpGetErrorCode(g_Handle);
			mexPrintf("Error getting Signal Format, error code = %d", ErrorCode);
			return;
		}
		
		g_NumberOfChannels = psf->Elements;
		g_BytesPerSample = g_NumberOfChannels * sizeof(long);
		
		g_SampleRateInMilliHz = 1000 * 1000;
		g_SignalBufferSizeInSamples = g_SampleRateInMilliHz/1000;
		g_SampleRateInHz = g_SampleRateInMilliHz/1000;
		
		if(fpSetSignalBuffer(g_Handle, &g_SampleRateInMilliHz,&g_SignalBufferSizeInSamples) != TRUE) 
		{
			mexPrintf("Error Setting Signal Buffer");
			return;
		}
		unsigned int SignalStrength, NrOfCRCErrors, NrOfSampleBlocks;
		
		Status = fpGetConnectionProperties(g_Handle, &SignalStrength,&NrOfCRCErrors, &NrOfSampleBlocks);
		if(!Status) {
			mexPrintf("Error acquiring Connection Properties");
			return;
		}
		
		g_startTime =  getUnixTimeStamp(); 
		g_CurCount = g_startTime;
		fpSetRefCalculation(g_Handle,1);
		if(!fpStart(g_Handle))
		{
			mexPrintf("Error starting recording");
			return;
		}
		DWORD dummy;
					  hObtainThread = CreateThread( 
						NULL,                   // default security attributes
						0,                      // use default stack size  
						threadObtain,       // thread function name
						NULL,          // argument to thread function 
						0,                      // use default creation flags 
						&dummy);   // returns the thread identifier 
		g_bIsConnected = true;				
}



	static SOCKET markerPassiveSocket;
	static SOCKET markerActiveSocket;
DWORD WINAPI threadMarkerServer( LPVOID lpParam ) 
{ 
	char curMarker[256];
	FILETIME* ft = new FILETIME;
	SYSTEMTIME* st = new SYSTEMTIME;
			
		markerPassiveSocket = socket(AF_INET,		
									SOCK_DGRAM,   	
									0);		
		if(markerPassiveSocket==INVALID_SOCKET) {
			 mexErrMsgTxt("bbci_acquire_en: Init. Could not create Passive socket\n");
			return 1;
		}
		SOCKADDR_IN serverInfo;
		serverInfo.sin_family = AF_INET;
		serverInfo.sin_addr.s_addr = INADDR_ANY;	
		serverInfo.sin_port = htons(1206);		
		
		if(bind(markerPassiveSocket, (LPSOCKADDR)&serverInfo, sizeof(struct sockaddr))==SOCKET_ERROR) {
			mexErrMsgTxt("Could not bind socket. Make sure that the port is free");
			return 1;
		}
		struct sockaddr_in si_other;
		int slen=sizeof(si_other);
		int rBytes;
		while(1/*(rBytes=recv(markerActiveSocket,curMarker, 256,0))!=SOCKET_ERROR*/) {
        rBytes=recvfrom(markerPassiveSocket,curMarker, 256,0, (struct sockaddr *) &si_other, &slen);
			if(rBytes!=SOCKET_ERROR)
			{
			
			if(!strncmp(curMarker,"QUIT_CALLED",rBytes)) {
				mexPrintf("QUITTING");
				break;
			}
			curMarker[4]=0;
			WaitForSingleObject(ghMutexMarkers, INFINITE );
			string curString = curMarker;
			mexPrintf("Received: %s\n",curMarker);
			markerData mData;
			mData.description = curString; 
			// mData.timeStamp = (unsigned long long) getUnixTimeStamp();
			//unsigned long long* recvTS = ((unsigned long long*) (&curMarker[248])); 
			//mexPrintf("%lld",*recvTS);
			mData.timeStamp = getUnixTimeStamp(); 
			bool isDone=false;
			unsigned long long lastTime = mData.timeStamp;
			
			unsigned long long startTime = getUnixTimeStamp();
			while(!isDone)
			{
				WaitForSingleObject(ghMutexCount, INFINITE );
					
				if(lastTime<=(g_CurCount))
				{
					isDone=true;
					mData.timeStamp=g_CurCount;
					//mexPrintf("\nlastTime:%lld CurCount:%lld\n",lastTime,g_CurCount);
				}
				ReleaseMutex( ghMutexCount);
			}
			unsigned long long endTime = getUnixTimeStamp();
			unsigned long long timeNeeded = endTime - startTime;
			
			gMarkerQueue.push(mData);
			ReleaseMutex(ghMutexMarkers);
			}
		
		//closesocket(markerActiveSocket);
		closesocket(markerPassiveSocket);
	}
	return 0;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxChar *pi;
		
	 unsigned char mac[6] = {0xcd, 0x02, 0x4c, 0x80, 0x07, 0x00};
	// Function case 1. 
	if(nrhs==2&&mxIsChar(prhs[0]) && mxIsStruct(prhs[1]))	{
		int len = mxGetM(prhs[0]) * mxGetN(prhs[0])+1;
		char* cPi = (char*) mxCalloc(len,sizeof(char));
		
		mxGetString(prhs[0],cPi,len);
		

		if(nlhs!=1)
		{
			mexErrMsgTxt("bbci_acquire_en: invalid number of outputs. init has only one output variable");
		}
		else if(!strcmp(cPi,"init")) 
		{
			if(g_bIsConnected) {
				mexPrintf("ALREADY CONNECTED. CLOSING PREVIOUS CONNECTION");
				quitTMSI();
			}
			
				gMarkerQueue = queue<markerData>();
				gDataQueue = queue<chData>();	
				mxArray* numChannels = mxGetField(prhs[1], 0,FIELD_Channels);
				if(numChannels)
					{
						double* t= (double*)mxGetData(numChannels);
						g_numCh = *t;
					}
				ghMutexData = CreateMutex(NULL,              // default security attributes
										  FALSE,             // initially not owned
										  NULL);             // unnamed mutex
				ghMutexMarkers = CreateMutex(NULL,              // default security attributes
										  FALSE,             // initially not owned
										  NULL);             // unnamed mutex
				ghMutexCount= CreateMutex(NULL,              // default security attributes
										  FALSE,             // initially not owned
										  NULL);             // unnamed mutex
		 
				mxArray* OUT_STATE;
				OUT_STATE = mxDuplicateArray(prhs[1]);
				plhs[0] = OUT_STATE;
				DWORD dummy;
						  hServerThread = CreateThread( 
							NULL,                   // default security attributes
							0,                      // use default stack size  
							threadMarkerServer,       // thread function name
							NULL,          // argument to thread function 
							0,                      // use default creation flags 
							&dummy);   // returns the thread identifier 
							
				
				initTMSI();
			
		}
		else if((!strcmp(cPi,"quit"))||(!strcmp(cPi,"close"))) {	
			quitTMSI();
		}
		else {
			mexErrMsgTxt("bbci_acquire_en: invalid string in first parameter. Did you mean init?");
		}
	}
	
	// Function case 2
	else if(nrhs==1&&nlhs==1 && mxIsStruct(prhs[0])) {
			WaitForSingleObject(ghMutexData, INFINITE );
			WaitForSingleObject(ghMutexMarkers, INFINITE );
			int count = gDataQueue.size();
			int count2 = gMarkerQueue.size();
			double* output = new double[count*g_numCh];
			//mexmexPrintf("COUNT: %d",count);
			for(int i=0;i<count;i++) {
				chData newData = gDataQueue.front();
				for(int j=0;j<g_numCh;j++)
				{
					output[i*g_numCh+j] = newData.channels[j];
				}
				gDataQueue.pop();
			}
			
			g_CurCount=0;
			gMarkerQueue = queue<markerData>();
			
			ReleaseMutex( ghMutexData);
			ReleaseMutex( ghMutexMarkers);
			plhs[0] = mxCreateDoubleMatrix(g_numCh,count,mxREAL);
			double  *start_of_output;
			start_of_output = (double *)mxGetPr(plhs[0]);
			memcpy(start_of_output, output, g_numCh*sizeof(double)*count);
			
			delete[] output;
			
		
	}
	// function case 3 and 4
		else if(nrhs==1&&nlhs>1 && mxIsStruct(prhs[0])) {

			bool isDone = false;
			WaitForSingleObject(ghMutexData, INFINITE );
			int preCount=gDataQueue.size();
			ReleaseMutex( ghMutexData);
			
			if(preCount>0)
			{
				WaitForSingleObject(ghMutexMarkers, INFINITE );
				WaitForSingleObject(ghMutexData, INFINITE );
				int count_markers = gMarkerQueue.size();
				double* output2 = new double[count_markers];
				char** output3 = new char*[count_markers];
				double* output4 = new double[count_markers];
				
				for(int i=0;i<count_markers;i++)
					output3[i] = new char[256];
				markerData* allMarkers =  new markerData[count_markers]; 
				for(int i=0;i<count_markers;i++)
				{
					allMarkers[i] = gMarkerQueue.front();
					strcpy(output3[i],gMarkerQueue.front().description.c_str());
					gMarkerQueue.pop();
					//mexmexPrintf("%s\n",output3[i]);
				}

				int count = gDataQueue.size();
				
				double* output = new double[count*g_numCh];

				
				unsigned long long* dataTime = new unsigned long long[count];
				

				
				for(int i=0;i<count;i++) {
					chData newData = gDataQueue.front();
					
					for(int j=0;j<g_numCh;j++)
					{
//						if(!newData.channels[j])
//							mexPrintf("\nZERO CH: %d", j);
						//output[j*count+i] = (newData.channels[j]/1000.0l); 
						output[j*count+i] = (newData.channels[j]); 
					
					}
					dataTime[i] = (unsigned long long)newData.timeStamp;
					gDataQueue.pop();
				}
				
				
				
				unsigned long long startTime = dataTime[0];
				for(int i=0;i<count_markers;i++)
				{
					// mexmexPrintf("\nStart time:%lld\n",startTime);
					unsigned long long diff = allMarkers[i].timeStamp - startTime;
					
//					fprintf(g_Fp, "%lld\n",diff);
					if(startTime>allMarkers[i].timeStamp)
					{	
					diff = 0;
					}

					output2[i] = diff;
				}
				
				plhs[0] = mxCreateDoubleMatrix(count,g_numCh,mxREAL);
				double  *start_of_output;
				start_of_output = (double *)mxGetPr(plhs[0]);
				memcpy(start_of_output, output, g_numCh*sizeof(double)*count);
				
				plhs[1] = mxCreateDoubleMatrix(1,count_markers,mxREAL);
				
				start_of_output = (double *)mxGetPr(plhs[1]);
				memcpy(start_of_output, output2, sizeof(double)*count_markers);
						
				int dims[2];
				dims[0] = 1;
				dims[1] = count_markers;
				
				//plhs[2] = mxCreateCellArray (2,dims);
				plhs[2] = mxCreateDoubleMatrix(1,count_markers,mxREAL);
				
				plhs[3] = mxDuplicateArray(prhs[0]);
				
				/*
				for(int i=0;i<count_markers;i++)
				{
				int arr[2];
				arr[0] = 0;
				arr[1] = i;
				int index = mxCalcSingleSubscript(plhs[2] , 2, arr);
				mxArray* curString =  mxCreateString(output3[i]);
				mxSetCell(plhs[2], index, curString);
				}
				*/
				for(int i=0;i<count_markers;i++)
				{
					string s(&output3[i][1]);
					istringstream(s) >> output4[i];
				}
				start_of_output = (double *)mxGetPr(plhs[2]);
				memcpy(start_of_output, output4, sizeof(double)*count_markers);
				
				delete[] allMarkers;
				delete[] output;
				delete[] output2;
				delete[] output4;
				delete[] dataTime;
				for(int i=0;i<count_markers;i++)
					delete[] output3[i];
				delete[] output3;
				ReleaseMutex( ghMutexMarkers);
				ReleaseMutex( ghMutexData);
				Sleep(3);
			}
			else
			{
				plhs[0] = mxCreateDoubleMatrix(g_numCh,0,mxREAL);
				plhs[1] = mxCreateDoubleMatrix(g_numCh,0,mxREAL);
				int dims[2];
				dims[0] = 1;
				dims[1] = 0;
				
				plhs[2] = mxCreateCellArray (2,dims);
				
				
				plhs[3] = mxDuplicateArray(prhs[0]);
			}
		
	
	}
	// Close
	else if(nlhs==0&&nrhs==1)
	{
		quitTMSI();
	}
}


void quitTMSI() {
	if(g_bIsConnected) {
		g_bIsConnected = false;
		g_numCh = 8;
		
		//TerminateThread(hServerThread,0);
		sockaddr_in si_other;
		SOCKET s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
		memset((char*) &si_other,0,sizeof(si_other));
		si_other.sin_family = AF_INET;
		si_other.sin_port = htons(1206);
		si_other.sin_addr.S_un.S_addr = inet_addr("127.0.0.1");
		char quitmsg[] = "QUIT_CALLED";
		int slen = sizeof(si_other);
		sendto(s,quitmsg, strlen(quitmsg),0, (sockaddr*) &si_other, slen);
		
		CloseHandle(hServerThread);
			//TerminateThread(hObtainThread,0);
		bTerminate=true;
		Sleep(1000);
		CloseHandle(hObtainThread);
		closesocket(markerPassiveSocket);
				
			
		WaitForSingleObject(ghMutexData, INFINITE );
		WaitForSingleObject(ghMutexMarkers, INFINITE );
		gDataQueue = queue<chData>();
		gMarkerQueue = queue<markerData>();
		ReleaseMutex( ghMutexMarkers);
		ReleaseMutex( ghMutexData);
		
		if(g_Handle) {
			fpStop(g_Handle);
			fpClose(g_Handle);
			Sleep(1000);
		}
		
		if( DeviceList != NULL ) 
			fpFreeDeviceList( g_Handle, NrOfDevices, DeviceList );
			
		DeviceList = NULL;
		if(g_Handle) {
			fpLibraryExit( g_Handle );
			g_Handle = NULL;
		}
			
		if(g_LibHandle) {
			FreeLibrary(g_LibHandle);
			g_LibHandle = NULL;
		}
	}
}