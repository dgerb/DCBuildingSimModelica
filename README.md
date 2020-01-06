# DCBuildingSimModelica
Python/Modelica platform for parametrically simulating DC buildings using NREL PVWatts solar profiles and PNNL Reference building load profiles


Requires Buildingspy and Modelica buildings library.
https://simulationresearch.lbl.gov/modelica/buildingspy/
https://simulationresearch.lbl.gov/modelica/

You might have to update the Modelica Buildings library reference to the current version.  To do that, open the package DC_Tool, and edit the line:
annotation (uses(Modelica(version="3.2.1"), Buildings(version="3.0.0")));
