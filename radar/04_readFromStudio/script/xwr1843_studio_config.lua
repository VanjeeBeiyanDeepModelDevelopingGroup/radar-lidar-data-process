-- Connection
-- set
ar1.FullReset()
ar1.SOPControl(2)
-- connect
ar1.Connect(9,921600,1000)
ar1.Calling_IsConnected()
ar1.SelectChipVersion("AR1642")
ar1.SelectChipVersion("XWR1843")
--ar1.SaveSettings('C:\\Users\\Administrator\\AppData\\Roaming\\RSTD\\ar1gui.ini')
-- load firmware
ar1.DownloadBSSFw("C:\\ti\\mmwave_studio_02_01_01_00\\mmWaveStudio\\Scripts\\..\\..\\rf_eval_firmware\\radarss\\xwr18xx_radarss.bin")
ar1.GetBSSFwVersion()
ar1.GetBSSPatchFwVersion()

ar1.DownloadMSSFw("C:\\ti\\mmwave_studio_02_01_01_00\\mmWaveStudio\\Scripts\\..\\..\\rf_eval_firmware\\masterss\\xwr18xx_masterss.bin")
ar1.GetMSSFwVersion()
-- SPI
ar1.PowerOn(0, 1000, 0, 0)
-- RF Power up
ar1.RfEnable()
ar1.GetMSSFwVersion()
ar1.GetBSSFwVersion()
ar1.GetBSSPatchFwVersion()

-- StaticConfig
-- Basic Configuration
ar1.ChanNAdcConfig(1, 1, 1, 1, 1, 1, 1, 2, 1, 0)
-- Advanced Configuration
ar1.LPModConfig(0, 0)
-- RF Init
ar1.RfInit()

-- Data Config
ar1.DataPathConfig(513, 1216644097, 0)
ar1.LvdsClkConfig(1, 1)
ar1.LVDSLaneConfig(0, 1, 1, 0, 0, 1, 0, 0)

-- Sensor Config
ar1.ProfileConfig(0, 76.5, 5, 4.8, 60, 0, 0, 0, 0, 0, 0, 4.007, 0, 256, 4652, 0, 0, 48)
ar1.ChirpConfig(0, 0, 0, 0, 0, 0, 0, 1, 1, 1)
ar1.FrameConfig(0, 0, 8, 128, 40, 0, 0, 1)

-- Set up DCA1000
ar1.GetCaptureCardDllVersion()
ar1.SelectCaptureDevice("DCA1000")

ar1.CaptureCardConfig_EthInit("192.168.33.30", "192.168.33.180", "12:34:56:78:90:12", 4096, 4098)
ar1.CaptureCardConfig_Mode(1, 2, 1, 2, 3, 30)
ar1.CaptureCardConfig_PacketDelay(25)
ar1.GetCaptureCardFPGAVersion()

-- Capture and Post Processing
ar1.CaptureCardConfig_StartRecord("C:\\ti\\mmwave_studio_02_01_01_00\\mmWaveStudio\\PostProc\\adc_data.bin", 1)
RSTD.Sleep(1000)

ar1.StartFrame()
RSTD.Sleep(5000)

ar1.StartMatlabPostProc("C:\\ti\\mmwave_studio_02_01_01_00\\mmWaveStudio\\PostProc\\adc_data.bin")
RSTD.Sleep(10000)