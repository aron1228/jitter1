#ifndef TMUJITTER_H_
#define TMUJITTER_H_

#include "TMUTask.h"
#include <float.h>

class TmuJitterUtil
{
public:

  /*
   *----------------------------------------------------------------------*
   * Struct: JitterParameter
   *
   * Purpose: Container to store jitter measurement parameters
   *
   *----------------------------------------------------------------------*
   * Description:
   *   The members in this structure can be categorized into 3 groups
   *     1. specified by a user
   *     2. specified by a user, possibly modified by processParameters()
   *     3. set by processParameters()
   * 
   *   The following parameters belong to group 1.
   *
   *   STRING pins:                 {@ | pin and/or pin group list}
   *     Name of pin(s) and/or pin group(s) to be measured
   *     Valid pins: all digital pins.
   *   DOUBLE datarate                     {}
   *     # max expected datarate of signal to be maesured
   *     # if omitted prescaler and deglitcher are used
   *   INT     prescaler;            {0..16}
   *   INT samples_per_meas:             {}
   *     # of samples to be taken
   *   INT ftstResult:               {0 | 1}
   *     Flag to specify whether taking functional test result into account
   *     for pass/fail judging and datalogging.
   *   
   *   The following parameters belong to group 2.
   *
   *   The following parameters belong to group 3.
   *
   *   TMU::EDGE_TYPE slope:
   *     Event slope of event . It's determined from the mode parameter
   *
   * Note:
   *   When you create an instance of this class, appropriate default values
   *   except pins are set to each variables in the class.
   *
   *----------------------------------------------------------------------*
   */
  struct JitterParameter
  {
      // Parameters specified by a user
      STRING  pins;
      STRING  mode;
      DOUBLE  datarate;
      DOUBLE  ui_ns;
      UINT32  prescaler;
      UINT64  samples_per_meas;
      Boolean exitOnNotFinished;
      DOUBLE  waitTimeout;
      Boolean ftstResult;
      Boolean histogram;

      // Parameters specified by a user, possibly modified by processParameters

      // Parameters set by processParameters()
      TMU::EdgeType slope;
      
      STRING testname;
      UINT32 allowed_outliers;
      DOUBLE outlier_value_limit;
      UINT32 outlier_loglevel;
      DOUBLE outlier_maxError;

      // Default constructor
      // all parameters are set to default values
      JitterParameter()
      : pins(""),
        mode(""),
        datarate(0.0),
        ui_ns(0.0),
        prescaler(0),
        samples_per_meas(1),
      exitOnNotFinished(FALSE),
        waitTimeout(0.0),
        ftstResult(FALSE),
        histogram(FALSE),
        slope(TMU::RISE),
        testname(""),
        allowed_outliers(0),
        outlier_value_limit(INFINITY),
        outlier_loglevel(0),
        outlier_maxError(0.0)
        {}
  };
  
  /*
   *----------------------------------------------------------------------*
   * Struct: JitterResult
   *
   * Purpose: Container to store jitter measurement results
   *
   *----------------------------------------------------------------------*
   */

  class JitterWaves   //储存jitter测量产生的结果
  {
  public:
    ARRAY_D wave_meas;
    Wave<DOUBLE> wave_jitter;
    Wave<INT>    wave_histo;
    std::vector < int > outlierindexes;
    std::vector < DOUBLE > outlier_error_values;
    DOUBLE min;
    DOUBLE max;
    DOUBLE stddev;
    DOUBLE rms;
    UINT32 outliers;        // number of detected outliers
    Boolean outlierfail;      // true if outliers > maxoutliers
  };

  class JitterResult : public TMU_TML_RESULT_BASE
  {
  public:
  Boolean periods;
    TMU::ApplicationType appType;
    Boolean histogram;
    UINT64 excpectedNumSamples;
    UINT32 maxoutliers;       // maximum allowed number of outliers
    UINT32 outlierloglevel;     // loglevel for outlier  0 = don't log: 1 = log occurrence: 2 = log detailed
    DOUBLE maxerrorvalue;     // maximum allowed sample error for outlier detection
    std::map < STRING, JitterWaves > pinWaveDatamap;
    JitterResult() : TMU_TML_RESULT_BASE(), periods(FALSE), appType(TMU::APP_RAW),histogram(FALSE),excpectedNumSamples(0),maxoutliers(0),outlierloglevel(0),maxerrorvalue(0.0) {}
    ~JitterResult()
    {
      pinWaveDatamap.clear();
    }

  };
  

  enum MinMaxEtc{
    RMS,
    P2P,
    MAX,
    MIN,
    STDDEV
  };
  
/*
 *----------------------------------------------------------------------*
 * Utility functions for Jitter measurements
 *----------------------------------------------------------------------*
 */

  /*
   *----------------------------------------------------------------------*
   * Routine:TmuJitterUtil::processParameters
   *
   * Purpose: Store given measurement parameters to the specified
   *          placeholder and determine addtional parameters
   *          which is necessary to execute measurement.
   *          Also performs some error checks on parameters.z
   *
   *----------------------------------------------------------------------*
   * Description:
   *   const STRING& pins:
   *     Measurement mode.
   *   DOUBLE datarate                     {}
   *     # max expected datarate of signal to be maesured
   *     # if omitted prescaler and deglitcher are used
   *   INT     prescaler;            {0..16}
   *
   *   INT samples_per_meas:             {}
   *     # of samples to be taken
   *   DOUBLE  treshold_a;               {V}
   *   INT ftstResult:               {0 | 1}
   *     Those are parameters for jitter measurement.
   *     See the descriptions in JitterParameter definition.
   *   JitterParameter& params:
   *     Container to hold parameters for the measurement.
   *   
   * Note:
   *
   *----------------------------------------------------------------------*
   */

  static void processParameters(const STRING& pins,                     //校验用户输入的参数，填充到JitterParameter结构体，并计算衍生参数。
                                                                   //这里的参数列表是输出，与前面那段不同，前面那段是用户输入的参数， 这里的是用来输出经过process运算过后的结果
                  const STRING& mode,
                                DOUBLE datarate,
                                DOUBLE UI_ns,
                                UINT32 prescaler,
                                UINT64 samples_per_meas,
                                INT exitOnNotFinished,
                                DOUBLE waitTimeout,
                                INT ftstResult,
                                Boolean histogram,
                                const STRING& testname,
                                UINT32 allowed_outliers,
                                UINT32 loglevelOutlier,
                                DOUBLE maxSampleError,
                                JitterParameter& params)
  {
    // If no pin for pins is specified, an exception is thrown.
     if (pins.size() == 0)
     {
       throw Error("TmuJitterUtil::processParameters()",
                     "Empty pins parameter.");
     }

     if ((datarate <= 0.0) && (prescaler == 0))
     {
       throw Error("TmuJitterUtil::processParameters()",
                    "Wrong datarate and prescaler parameter.");
     }

     if (((mode == "DATA") || (mode == "DATA_PERIOD"))  && (UI_ns <= 0.0))
     {
       throw Error("TmuJitterUtil::processParameters()",
                   "valid UI value expected for jitter mode DATA");
     }


     // TODO replace with tmu own definitions
       params.ui_ns = UI_ns;
       params.mode  = mode;
       params.slope = TMU::RISE;           //设置事件触发条件，jitter测试应设置成rise——fall
       params.pins  = pins;
       params.datarate = datarate;
       params.samples_per_meas = samples_per_meas;
       params.ftstResult = (ftstResult == 0)? FALSE:TRUE;    //标志，指定是否将功能测试结果纳入账户用于合格/不合格判断和数据记录。
       params.prescaler = prescaler;
       params.histogram = histogram;
       params.exitOnNotFinished = (exitOnNotFinished == 1)?TRUE:FALSE;
       params.waitTimeout   = waitTimeout;
       params.testname = testname;
       params.allowed_outliers = allowed_outliers;
       params.outlier_loglevel = loglevelOutlier;
       /*异常值处理：
           如果用户提供了 maxSampleError，就直接使用它作为异常值限制。
           如果用户未提供，就从测试限制表中加载限制。
           调整允许的异常值数量：根据样本数和用户输入的比例计算允许的异常值数量。 */
       if(maxSampleError > -1) // if maxSampleError > -1 second use it 
       {
         params.outlier_value_limit = maxSampleError;
       }
       else
       {
         // load the limit for p2pjitter to be able to detect outliers in calulateJittervalues
         string testsuiteName,errormessage,JitterTestnameP2P("p2pJitter");
         vector<string> names;
         if(extractTestNameVector(params.testname,names,5))
         {
           JitterTestnameP2P  =  names.at(1);
         }
         Boolean sIsLimitTableLoaded = V93kLimits::tmLimits.loadOnDemandAndReturnStatus();

         LIMIT limit;
         GET_TESTSUITE_NAME(testsuiteName);
         DOUBLE factor = 1.0;
         getTmuLimits(factor,sIsLimitTableLoaded,
             testsuiteName,JitterTestnameP2P,"ns","1[s]",errormessage,limit);
         DOUBLE dlimit;
         limit.getHigh(&dlimit);
         params.outlier_value_limit = dlimit;
       }

       if(params.allowed_outliers > 0)
       {
         params.allowed_outliers = (UINT32) ceil (params.samples_per_meas / 1000.0 * params.allowed_outliers);//params.allowed_outliers用户输入的是百分比，
                                                                                                              // 然后通过这个函数处理后再转伟32位无符号整数储存
       }

  }
  
  /*
   *----------------------------------------------------------------------*
   * Routine: TmuJitterUtil::doMeasurement
   *
   * Purpose: Perform setup, execution and result retrieval
   *          for rise/fall time measurement with per pin TIA
   *
   *----------------------------------------------------------------------*
   *
   * Description:
   *   const JitterParameter& params:
   *     Container to hold parameters for the measurement.
   *   JitterResult& results:
   *     Container to store measurement results
   *
   * Note:
   *   'static' variables are used in this function to keep some information
   *   which is refered in the execution for sites where ON_FIRST_INVOCATION 
   *   block is not executed.
   *
   *----------------------------------------------------------------------*
   */
  static void doMeasurement(const JitterParameter& params,
                            JitterResult& results)
  {
      // Measurement setup, execution and result retrieval are done
      // through task object
      TMU_TASK task;
      TMU_RESULTS siteResults;
      results.periods = FALSE;
      if((params.mode == "PERIOD")
      || (params.mode == "DATA_PERIOD"))
      {
        results.periods = TRUE;
      }
      INT site;
      STRING ErrorString;
      Boolean exception_caught = FALSE;
      Boolean finished = TRUE;
      results.histogram = params.histogram;
      results.maxoutliers = params.allowed_outliers;
      results.outlierloglevel = params.outlier_loglevel;
      results.maxerrorvalue = params.outlier_value_limit / 2.0;


      try
      {
        CONNECT();
        site = CURRENT_SITE_NUMBER();
        ON_FIRST_INVOCATION_BEGIN(); //多站点测试时，可以并行测试，单site测试就忽略
        results.errorStr.clear();
        //avoid ON_FIRST_INVOCATION_END executed so catch exceptions
        try
        {
          if(params.datarate > 0.0)
          {
            task.pin(params.pins).setDatarate(params.datarate MBps);     //如果输入的是datarate，setDatarate设置数据速率，将自动计算预分频和丢弃的值
          }
          else if((params.prescaler > 0))
          {
            task.pin(params.pins)
                          .setPreScaler(params.prescaler)
                          .setInterSampleDiscard(0);   //设置预分频数和初始化丢弃数
          }


          task.pin(params.pins).setEdgeSelect(params.slope)   //设置触发事件的边缘类型。
                                   .setNumSamples(params.samples_per_meas)//每次测量的采样数
                                   .setInitialDiscard(1)//设置丢弃数
                                   .setNumMeasurements(1)
                                   .setNumShots(1);

          /*
           * TMU measurement execution
           */
          task.setAsPrimary();
          task.execute();

          task.pin(params.pins).waitTMUfinished(params.waitTimeout,finished);

          if((params.exitOnNotFinished == FALSE)
              || (finished == TRUE))
          {
            task.pin(params.pins).uploadRawResults(siteResults);//上传原始数据
          }

           //异常处理及报告
        }
        catch(Error& e)  // catch here and set errorflag so     ON_FIRST_INVOCATION_END will be called
        {
          ErrorString = (e.msg()).c_str();
          exception_caught = TRUE;
        }

        ON_FIRST_INVOCATION_END();
        results.funcTestResult  = task.getStatus();      //•返回上次执行的状态。它可以是SUCCESS， FAILED， ABORTED或not_performed。
        results.pinlist = params.pins;
        results.excpectedNumSamples = params.samples_per_meas;
        /*
         * TMU measurement result retrieval
         */
        if((params.exitOnNotFinished == TRUE)
            && (finished == FALSE))
        {
          throw TESTMETHOD_API::TMException("TmuJitterUtil::doMeasurement()","TMU has not finished","");
        }
        if(exception_caught)
        {
          throw TESTMETHOD_API::TMException("TmuJitterUtil::doMeasurement()",ErrorString.c_str(),"");
        }

      }
      catch(TESTMETHOD_API::TMException& e)
      {
        STRING error = e.getOrigin() + " : ";
        error += e.msg().c_str();
        TESTSET().cont(TM::CONTINUE).
            judgeAndLog_ParametricTest("", "", TM::Fail, 0);
        throw Error("TmuJitterUtil::doMeasurement() ",error);
      }

      if (FALSE == siteResults.checkSiteAppTypeStored(TMU::APP_RAW))      //检查是否有测量结果，没有就报错
      {
       results.noResultsForSite = TRUE;
       char buf[128];
       sprintf(buf,"No jitter measurement results  available for site %d",site);
       results.errorStr += buf;
      }
      else
      {
        // Get the status result from functional test if specified.
        results.funcTestResult = task.getStatus();//获取最后执行的状态

        STRING_VECTOR& pinvector = results.pinvector;
        if (siteResults.getNumberOfMeasuredPins(pinvector, results.appType,site) == 0)
        {
          results.noPinsHasResults = TRUE;
          char buf[128];
          sprintf(buf,"No pinresults available for Jitter measurement for site %d",CURRENT_SITE_NUMBER() );
          results.errorStr += buf;
          return;
        }
        getMissingPins(results.missingpins,pinvector,results.pinlist);

        for(size_t pinindex = 0; pinindex < pinvector.size(); pinindex++)
        {

          STRING& pin = pinvector[pinindex];
          UINT32 maes = 1; // armingID
          JitterWaves& waves = results.pinWaveDatamap[pin];
          siteResults.getPinMeasuredResults(pin,waves.wave_meas,TMU::APP_RAW,maes);
          if(calculateJitterValues(waves.wave_meas,waves,params))    //calculateJitterValues(ARRAY_D& wave_meas,JitterWaves& wave,const JitterParameter& params)
          {
            INT size = waves.wave_jitter.getSize();

            waves.min    = waves.wave_jitter.min();
            waves.max    = waves.wave_jitter.max();
            waves.stddev = waves.wave_jitter.stddev();
            waves.rms  = waves.wave_jitter.rms();

            DOUBLE jittermean = waves.wave_jitter.mean();  // determine zero point
            for(int ii = 0; ii < size;ii++)
            {
              waves.wave_jitter[ii] -= jittermean;         // allign around zero
            }
            waves.min -= jittermean;                   // allign around zero
            waves.max -= jittermean;
            if(results.histogram == TRUE)
            {
              INT nBins = (INT) (floor(sqrt((DOUBLE)(size - 1))) - 1);
              waves.wave_histo.histogram(waves.wave_jitter,nBins ,waves.min,waves.max);
            }
          }
        }
      }
  }

  /*
   *----------------------------------------------------------------------*
   * Routine: TmuJitterUtil::judgeAndDatalog
   *
   * Purpose: Perform pass / fail judgement and datalogging
   *          for jitter measurement
   *
   *----------------------------------------------------------------------*
   * Description:
   *   JitterResult& results:
   *     Container to store measurement results
   *   INT jitterHistogram:          {0 | 1}
   *     Flag to specify whether jitter histogram is logged
   *   INT ftstResult:               {0 | 1}
   *     Flag to specify whether taking functional test result into account
   *     for pass/fail judging and datalogging.
   *   
   * Note:
   *   Jitter histogram is logged first.
   *   judgeAndLog_ParametricTest() is called for results in the order of
   *   TM::STDEV, TM::MINMAX, and finally functional test
   *   results. Even if one of results is FAIL, test method is executed
   *   until all avaiable results are logged.
   *
   *----------------------------------------------------------------------*
   */
  static void judgeAndDatalog(const string& testname, JitterResult& results,
                              INT ftstResult)
  {
    TMU_RESULTS siteResults;

    static bool sIsLimitTableLoaded = false;
    static bool sIsLimitTableUsed[5] = {false,false,false,false,false};
    static string JitterTestname[5] ={"","","","",""};
    static DOUBLE factor[5] = {1.0,1.0,1.0,1.0,1.0};
    static LIMIT limit[5];

    static bool   hasfailed = false;
    static string ErrorStr;
    string testsuiteName;
    GET_TESTSUITE_NAME(testsuiteName);

      ON_FIRST_INVOCATION_BEGIN();
      string errormessage;
      bool limitLoaded[5] = {false,false,false,false,false};
      hasfailed = false;
      ErrorStr.clear();
      //check whether limit table is used.
      vector<string> names;
      if(extractTestNameVector(testname,names,5))
      {
        JitterTestname[RMS]   = names.at(0);
        JitterTestname[P2P]   = names.at(1);
        JitterTestname[MAX]   = names.at(2);
        JitterTestname[MIN]   = names.at(3);
        JitterTestname[STDDEV]= names.at(4);
      }

      sIsLimitTableLoaded = V93kLimits::tmLimits.loadOnDemandAndReturnStatus();

      for (int i = RMS;i <= STDDEV; i++)
      {
        sIsLimitTableUsed[i] = sIsLimitTableLoaded;
        limitLoaded[i] = getTmuLimits(factor[i],sIsLimitTableUsed[i],
            testsuiteName,JitterTestname[i],"ns","1[s]",errormessage,limit[i]);
      }

      if (!limitLoaded[RMS] || !limitLoaded[P2P] || !limitLoaded[MAX]
          || !limitLoaded[MIN] || !limitLoaded[STDDEV])
      {
        hasfailed = true;
        ErrorStr += "\n Testsuite ";
        ErrorStr += testsuiteName + ":\n";
        ErrorStr += errormessage + "\n";
      }

    ON_FIRST_INVOCATION_END();
    if(hasfailed)
    {
      results.limitLoadError = true;
      results.errorStr += ErrorStr;
    }
    if(results.hasFailed())
    {
      TESTSET().cont(TM::CONTINUE).
          judgeAndLog_ParametricTest("", "", TM::Fail, 0);
      PUT_DATALOG(results.errorStr);
  		return;
    }
    try
    {
      logMissingPins(results.missingpins);
      STRING_VECTOR& pinvector = results.pinvector;
      for(size_t pinindex = 0; pinindex < pinvector.size(); pinindex++)
      {
        STRING& pin = pinvector[pinindex];
        Boolean failed = FALSE;
        JitterWaves& waves = results.pinWaveDatamap[pin];
        INT size = waves.wave_jitter.getSize();
        if(size <= 0)
        {
          failed = TRUE;
        }
        else
        {
          DOUBLE rms = waves.wave_jitter.rms();
          DOUBLE peek2peek = waves.max - waves.min;
          DOUBLE max = waves.max;
          DOUBLE min = waves.min;
          DOUBLE stddev = waves.stddev;
          LIMIT limitoutlier(TM::GE,0.0,TM::LE,(DOUBLE)results.maxoutliers);
          LIMIT limitoutlierval(TM::GE,-results.maxerrorvalue,TM::LE,results.maxerrorvalue);
          limitoutlierval.unit("ps;-12;s;0");
          if(results.outlierloglevel > 0)
          {
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,"Number of outliers",limitoutlier,waves.outliers);
            if((results.outlierloglevel == 2) && (!waves.outlierindexes.empty()))
            {
              char buf[256];
              memset(buf,'\0',256);
              sprintf(buf,"\n%s: outlier max error: %f to %f [ps]\n",pin.c_str(),-results.maxerrorvalue *1e12,results.maxerrorvalue *1e12);
              STRING logstr(buf);
              for(unsigned int ii = 0; ii < waves.outlierindexes.size();ii++)
              {
                sprintf(buf,"%s: outlier error at index %d = %f [ps]\n",pin.c_str(),waves.outlierindexes[ii],waves.outlier_error_values[ii]*1e12);
                logstr += buf;
              }
              PUT_DATALOG(logstr );
            }
          }
          if(sIsLimitTableUsed[RMS]) {// limit from testtable
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[RMS], V93kLimits::tmLimits,rms* factor[RMS]);
          } else {
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[RMS], limit[RMS],rms* factor[RMS]);
          }
          if(sIsLimitTableUsed[P2P]) {// limit from testtable
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[P2P], V93kLimits::tmLimits,peek2peek* factor[P2P]);
          } else {
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[P2P], limit[P2P],peek2peek* factor[P2P]);
          }
          if(sIsLimitTableUsed[MAX]) {// limit from testtable
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[MAX], V93kLimits::tmLimits,max* factor[MAX]);
          } else {
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[MAX], limit[MAX],max* factor[MAX]);
          }
          if(sIsLimitTableUsed[MIN]) {// limit from testtable
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[MIN], V93kLimits::tmLimits,min* factor[MIN]);
          } else {
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[MIN], limit[MIN],min* factor[MIN]);
          }
          if(sIsLimitTableUsed[STDDEV]) {// limit from testtable
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[STDDEV],V93kLimits::tmLimits,stddev* factor[STDDEV]);
          }
          else {
            TESTSET().cont(TM::CONTINUE).judgeAndLog_ParametricTest(pin,JitterTestname[STDDEV],limit[STDDEV],stddev* factor[STDDEV]);
          }
        }

        INT64 checksize = (results.periods == TRUE)?
            results.excpectedNumSamples -1: results.excpectedNumSamples;

        if(waves.wave_jitter.getSize() != checksize)
        {
          TESTSET().cont(TM::CONTINUE).reliability(TM::UNRELIABLE).
                judgeAndLog_ParametricTest(pin,testname, 0);
        }

      } // for(int pinindex ......

      // Finally judge and log the functional test result if available
      if (ftstResult)
      {
        if (results.funcTestResult == TMU::SUCCESS)
    	{
          TESTSET().cont(TM::CONTINUE).
              judgeAndLog_FunctionalTest(true);
    	}
    	else
    	{
    	  TESTSET().cont(TM::CONTINUE).
              judgeAndLog_FunctionalTest(false);
    	}
      }
    }
    catch(TESTMETHOD_API::TMException& e)
    {
      throw Error("TmuJitterUtil::judgeAndDatalog()",(e.msg()).c_str());
    }
  }

  /*
   *----------------------------------------------------------------------*
   * Routine: TmuJitterUtil::reportToUI
   *
   * Purpose: Output jitter measurement results to Report Window
   *
   *----------------------------------------------------------------------*
   * Description:
   *   JitterResult& results:
   *     Container to store measurement results
   *   INT ftstResult:               {0 | 1}
   *     Flag to specify whether to output functional test result
   *   const STRING& output          {None | ReportUI}
   *     Flag to specify whether to output 
   *     "None" means no output and "ReportUI" means to do output
   *
   * Note:
   *   "debug_analog" testflow flag should be set to see the jitter histogram
   *   on the singal analyzer
   *
   *----------------------------------------------------------------------*
   */
  static void reportToUI(JitterResult& results,
      const STRING& output,
      INT ftstResult)
  {
    TMU_RESULTS siteResults;
    try
    {
      // If output parameter is different from "ReportUI", just retun
      if (( output != "ReportUI" ) && ( output != "ReportUI_nosamples" )) {
        return;
      }
      printMissingPins(results.missingpins);

      // ------!!!----- allways avoid nullpointer access --------------
      // If nothing stored in results for this applicationtype return
      if(results.hasFailed())
      {
        printf(results.errorStr.c_str());
        return;
      }
      STRING_VECTOR& pinvector = results.pinvector;

      for(size_t pinindex = 0; pinindex < pinvector.size(); pinindex++)
      {
        STRING& pin = pinvector[pinindex];
        Boolean failed = FALSE;
        JitterWaves& waves = results.pinWaveDatamap[pin];
        INT size = waves.wave_jitter.getSize();
        if(size <= 0)
        {
          failed = TRUE;
          printf("Pin: [%s]  \n", pin.c_str());
          printf( "rmsJitter    **********   sec\t*** NO VALID RESULT ***\n");
          printf( "p2pJitter     **********   sec\t*** NO VALID RESULT ***\n");
        }
        else
        {
          printf("Pin: [%s]  \n", pin.c_str());
          if( output != "ReportUI_nosamples" )
          {
            for(int k = 0; k < size;k++)
            {
              printf("sampleJitter %-6d\t= % 13g   sec\n",k,waves.wave_jitter[k]);
            }
          }
          printf( "------------------------------------\n");
          printf( "rmsJitter          \t= % 13g   sec\n",waves.rms);
          printf( "p2pJitter          \t= % 13g   sec\n",waves.max - waves.min);
          printf( "maxJitter          \t= % 13g   sec\n",waves.max);
          printf( "minJitter          \t= % 13g   sec\n",waves.min);
          printf( "stddevJitter       \t= % 13g   sec\n",waves.stddev);
          printf( "------------------------------------\n");
          if(results.histogram == TRUE)
          {    // calculate bins and ARRAYs for signal analyzer
            INT nBins = (INT) (floor(sqrt((DOUBLE)(size - 1))) - 1);
            DOUBLE binWidth = (waves.max - waves.min) / nBins;
            ARRAY_D xseries;
            ARRAY_D histogram;
            xseries.resize(nBins);
            xseries[0] = waves.min;
            for(int i = 0;i < nBins;i++)
              xseries[i] = 1e12*(waves.min + (i*binWidth)); // bins as pico sec.

            histogram.resize(waves.wave_histo.getSize());
            for (int i = 0; i < waves.wave_histo.getSize(); i++)
              histogram[i] = 1.0*waves.wave_histo[i];
            PUT_DEBUG(pin.c_str(), "Jitter Histogram",xseries,histogram);

            Wave < DOUBLE > printwave(waves.wave_jitter);
            for (INT i = 0; i < waves.wave_jitter.getSize();i++)
              printwave[i] *= 1e12;
            PUT_DEBUG(pin.c_str(), "Edge Errors",printwave);
          }

        }
        if(results.outlierloglevel > 0)
        {
          if(waves.outlierfail)
          {
            printf( "Too many outliers detected: %d\n",waves.outliers);
          }
          if(results.outlierloglevel == 2)
          {
            printf( "max error = +-%f [ps]\n",results.maxerrorvalue * 1e12);
            for(unsigned int i = 0; i < waves.outlierindexes.size();i++)
            {
              printf( "sample error at index %d = %f [ps] \n",i,waves.outlier_error_values[i]* 1e12);
            }
          }
        }
        // are there all asumed samples measured
        INT64 checksize = (results.periods == TRUE)?
            results.excpectedNumSamples -1: results.excpectedNumSamples;
        if(size != checksize)
        {
          printf("           **********   not all samples measured\n");
        }

      } // for(int pinindex ......


      fflush(stdout);
    }
    catch(Error& e)
    {
      fflush(stdout);
      throw Error("TmuJitterUtil::reportToUI()",(e.msg()).c_str());
    }

  }
/*mode有四种，data:测量数据相对UI的抖动的，dataperiod：测量相邻信号间的抖动的，clk：测时钟的，period测相邻时钟抖动的





*/
  static Boolean calculateJitterValues(ARRAY_D& wave_meas,JitterWaves& wave,const JitterParameter& params)
  {
    Boolean ret = true;
    ARRAY_D& source = wave_meas;
    Wave<DOUBLE>& result = wave.wave_jitter;
    TMU_RESULTS siteResults;

    DOUBLE maxdiff = params.outlier_value_limit / 2;       //根据 params.outlier_value_limit 设置抖动值的上、下限阈值。
   // maxdiff 和 mindiff                                     //用于判断样本值是否为异常值。
    DOUBLE mindiff = -1.0 * maxdiff;
    wave.outliers = 0;

    TMU::DATARATE_DEPENDENT_VALUES_TYPE settings;   //这个是用来获得预分频和丢弃数
    TMU_TASK task;
    try
    {
      if(source.size() <= 1)                        //source就是wave.mea，里面存储的是原始波形数据
        throw TESTMETHOD_API::TMException("TmuJitterUtil","calculateJitterValues()",
            "source array empty or only one sample");

      DOUBLE dummy,modTmp,periodMean;
      periodMean = dummy = 0.0;
      // determine period
      INT size = source.size();
      if(size < 2) return false;
      DOUBLE halfperiod;
      task.pin(params.pins).getDatarateDependendSettings(settings);   //检索与数据相关的分频数和丢弃数
      if((params.mode == "DATA")
          || (params.mode == "DATA_PERIOD"))
      {
        periodMean = params.ui_ns * 2e-9;   //是单位时间间隔的两倍，并将单位换算成s，UI指的是一个单元间隔，不是指的一个周期间隔
        halfperiod = params.ui_ns * 1e-9;
      }
      else
      {
       	DSP_REG1(source,&periodMean,&dummy,0,size-1);  //最小二乘法，用来计算另一个模式的周期
        periodMean = periodMean / settings.prescaler / (DOUBLE)(settings.interdiscard+1);//tmu测量间隔是有最小值的，所以会通过将原高频信号先分频，再丢弃来延长时间，达到所需要
                                                                                        //最小值，所以实际周期就可以反过来
        halfperiod = periodMean / 2.0;
      }

      // compute sample differences periodMean
      if(params.mode == "DATA")
      {
        result.setSize(size );
        for(int i = 0;i < size ;i++)
        {
          modTmp = remainder(source[i]-source[0],halfperiod);     //计算当前信号在一个UI中的相对位置，remainder是取余函数
          result[i] = modTmp;
        }
      }
      else
      {
        if(params.mode == "CLK")
        {
          result.setSize(size );
          for(int i = 0;i < size ;i++)
          {
            modTmp = remainder(source[i]-dummy,periodMean);
            result[i] = (modTmp > halfperiod)? modTmp - periodMean : modTmp;//将wck的抖动控制在【-半周期，+半周期】

            if( (params.allowed_outliers > 0) && ((result[i] > maxdiff) || (result[i] < mindiff)) ) {
              wave.outliers++;
              wave.outlierindexes.push_back(i);
              if(params.outlier_loglevel == 2)
              {
                wave.outlier_error_values.push_back(result[i]);         /*每次检测到抖动值超出阈值时，wave.outliers 计数器递增。
                                                                       将异常值的索引（i）存储在 wave.outlierindexes 中。
                                                                       如果 params.outlier_loglevel 设置为 2，还将异常值的误差（result[i]）存储在 wave.outlier_error_values 中。
                */
              }
            }
          }
        }
        else if(params.mode == "PERIOD") //  -> PPTia correlating results    相邻信号的抖动
        {
          result.setSize(size -1);
          for(int i = 1;i < size ;i++)
          {
            modTmp = remainder(source[i]-source[i-1],periodMean);
            result[i-1] = modTmp;
            if( (params.allowed_outliers > 0) && ((result[i-1] > maxdiff) || (result[i-1] < mindiff)) ) {
              wave.outliers++;
              wave.outlierindexes.push_back(i);
              if(params.outlier_loglevel == 2)
              {
                wave.outlier_error_values.push_back(result[i]);
              }
            }
          }
        }
        else  // DATA_PERIOD -> PPTia correlating results            
        {
          result.setSize(size -1);
          for(int i = 1;i < size ;i++)
          {
            modTmp = remainder(source[i]-source[i-1],periodMean);
            result[i-1] = modTmp;
          }
        }
      }

      DOUBLE mean = result.mean();
      if(!wave.outlierindexes.empty() && (params.allowed_outliers > 0))
      {
        for( int idx = 0; idx < (int)wave.outlierindexes.size(); idx++ )
        {
          result[wave.outlierindexes[idx]] = mean; //吧异常的结果放到mean里去
        }
        if(wave.outliers > params.allowed_outliers)
        {
          wave.outlierfail = true;
        }
      }
    }
    catch(Error& e)
    {
      ret = FALSE;
    }
    catch(...)
    {
      ret = FALSE;
    }
    return ret;
  }


};
  
#endif /*TMUJITTER_H_*/