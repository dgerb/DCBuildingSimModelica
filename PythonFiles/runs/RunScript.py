# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
plt.close('all')
import os, sys
PROJECTPATH = os.path.dirname(os.path.dirname(__file__))
sys.path.append(PROJECTPATH)
import packageDC.DCBuildingsPython as dc
reload(dc)

def main():

    SIM = False
    PLOT = False
    SIMBRUCE = False
    SHOWGUI = False
    PLOTTRAN = False
    CLEAN = True

    SIM = True
    PLOT = True
    CLEAN = False
    
    if CLEAN:
        dc.cleanSimDirectories()
    
    # a = 'MOB_'
    # a = 'APT_'
    
    # b = 'SF_'
    # b = 'Miami_'
    # b = 'Duluth_'
    
    # c = 'Equip'
    # c = 'LED'
    # c = 'VFD'
    
    # runs = {a+b+c: ['Opportunity']}

    # a = '1'
    # a = '2'
    # a = '3'
    
    # b = 'AC'
    # b = 'DC'
    
    # runs = {'Test'+a+'_'+b: ['CERCTest']}

    # runs = {'CSUTest': ['CSUTest']}
    # runs = {'CSUSteve': ['CSUTest']}
    # runs = {'CSU_AC1': ['CSUTest']}
    # runs = {'CSU_AC2': ['CSUTest']}
    # runs = {'CSU_AC3': ['CSUTest']}
    # runs = {'CSU_DC': ['CSUTest']}
            
    # runs = {'Bruce_MOB_LA': ['Bruce']}
    
    # runs = {'MOB_SF': ['SteveMELs']}
    # runs = {'MOB_Miami': ['SteveMELs']}
    runs = {'MOB_Duluth': ['SteveMELs']}

    # runs = {'MOB_LA': ['SolarBattery', 'SolarBatteryFuture']}
    # runs = {'MOB_SF': ['SolarBattery', 'SolarBatteryFuture']}
    # runs = {'Rest_LA': ['SolarBattery', 'SolarBatteryFuture']}
    # runs = {'Rest_SF': ['SolarBattery', 'SolarBatteryFuture']}
    # runs = {'Retail_LA': ['SolarBattery', 'SolarBatteryFuture']}
    # runs = {'Retail_SF': ['SolarBattery', 'SolarBatteryFuture']}

    # runs = {'Hotel_LA': ['SolarBattery', 'SolarBatteryFuture']}
    # runs = {'Hotel_SF': ['SolarBattery', 'SolarBatteryFuture']}
    # runs = {'Apt_LA': ['SolarBattery']}
    # runs = {'Apt_SF': ['SolarBattery']}
    # runs = {'Market_LA': ['SolarBattery']}
    # runs = {'Market_SF': ['SolarBattery']}    
    
    if (not PLOT) and SIM:
        dc.cleanSimDirectories()
    for gridName in runs:
        dc.createProfileMat(gridName, PROJECTPATH)
        for experimentName in runs[gridName]:
            print '\n', gridName, experimentName
            mParams, pParams = dc.getSimParams(PROJECTPATH, gridName, \
                experimentName)
            converterList = dc.appendModelParams( \
                mParams, pParams, gridName, PROJECTPATH)
            simDirectory = dc.getSimDirectory(gridName, PROJECTPATH)
            # Simulating
            if SIM or SIMBRUCE:
                # calibrate before each simulation
                for runName in mParams:
                    # dc.bruceJson(mParams[runName], pParams[runName], \
                    #     PROJECTPATH, gridName, experimentName, runName)
                    runMParams = {runName: mParams[runName]}
                    runPParams = {runName: pParams[runName]}
                    # Calibrate
                    calibMParams, calibPParams = dc.getCalibrationParams( \
                        runName, runMParams, runPParams)
                    dc.appendModelParams( \
                        calibMParams, calibPParams, gridName, PROJECTPATH)
                    dc.appendPeakEfficiencies(calibMParams, gridName, \
                        simDirectory, converterList, PROJECTPATH)
                    dc.simulate(calibMParams, calibPParams, True, \
                        gridName, simDirectory, SHOWGUI)
                    # Run Sim
                    runMParams = dc.calibrate(gridName, converterList, \
                        runMParams, runPParams, PROJECTPATH, simDirectory)
                    if SIMBRUCE:
                        dc.bruceJson(runMParams[runName], runPParams[runName], \
                            PROJECTPATH, gridName, experimentName, runName, \
                            isCalibrated=True)
                    # dc.bruceJson(mParams[runName], pParams[runName], \
                    #     PROJECTPATH, gridName, experimentName, runName, \
                    #     isCalibrated=True)
                    if not SIM:
                        continue
                    dc.simulate(runMParams, runPParams, False, \
                        gridName, simDirectory, SHOWGUI)

            # Plotting
            # gridInfo is an array of [[gridNames], [gridColors]]
            gridInfo = [['AC', 'DC'], ['r', 'b']]
            if (PLOT):
                plotDirectory = dc.getPlotDirectory(gridName, PROJECTPATH)
                dc.plotTranAndStats(PLOTTRAN, experimentName, mParams, \
                    pParams, converterList, gridInfo, gridName, simDirectory, plotDirectory)
        dc.createProfileMat(gridName, PROJECTPATH)

    if CLEAN:
        dc.cleanSimDirectories()


if __name__ == '__main__':
    main()
