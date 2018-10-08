# OldSpot: A Pre-RTL Model for Fine-grained Aging and Lifetime Optimization
OldSpot is a simulation framework for estimating the lifetime of an SoC at the architectural unit level with support for failure tolerance. It can simulate heterogeneous systems with different types of cores (i.e. general-purpose cores or application specific accelerators) and resources shared among them such as last-level cache. For more detailed information on OldSpot, please refer to:
* A. Roelke, X. Guo, and M. R. Stan, "OldSpot: A Pre-RTL Model for Fine-grained Aging and Lifetime Optimization." In *Proceedings of the 36th IEEE International Conference on Computer Design*, October 2018.
* Chapter 3 of A. Roelke, "Improving Reliability and Security with Aging and Pre-RTL Modeling." PhD Thesis, University of Virginia, 2018.

OldSpot also incorporates device-level models which it uses to compute the aging of each microarchitectural units. It includes models for negative bias temperature instability (NBTI), electromigration (EM), hot-carrier injection (HCI), and time-dependent dielectric breakdown (TDDB).

If you use OldSpot in your research, we would appreciate a citation to:
```
@inproceedings{roelke2018oldspot,
  author = {A. Roelke and X. Guo and M. R. Stan},
  title = {{OldSpot}: A Pre-RTL Model for Fine-grained Aging and Lifetime Optimization},
  booktitle = {36th IEEE International Conference on Computer Design (ICCD)},
  year = {2018},
  month = {October},
  pages = {148--151}
}
```

## Compiling OldSpot
Simply go to the root of this project and type `make`.

### Dependencies
OldSpot requires the following dependencies:
* GCC 4.8+
* [tclap](http://tclap.sourceforge.net/) (Available on Debian systems as `libtclap-dev`)
* [pugixml](https://pugixml.org/) (Available on Debian as `libpugixml-dev`)

## Running OldSpot
After compiling, OldSpot can be executed by running `./oldspot [config file]`.

### Configuration File
The configuration file is an XML document describing failure propagation through the system. It contains two sections: one describing each architectural unit and one describing how failures propagate. In the first section, each unit's specification looks something like:
```
<unit type="unit" name="core">
    <default vdd="1.1" />
    <default temperature="350" />
    <default frequency="1000" />
    <redundancy type="serial" count="1" />
    <trace file="example/one.trace" failed="" />
</unit>
```
Each unit is specified using a `<unit>` element, which as two attributs:
* `type`: The type of unit being described. Valid values are `unit`, `core`, `logic`, and `memory`.
* `name`: A unique string identifier for the unit.
Each `<unit>` element can also have the following elements:
* `<default>`: Specifies a default value for a type of data if it is not present in the input trace, which is specified using an attribute. Each `<default>` tag only describes one type of data, but a `<unit>` element can have any number of `<default>` elements.  Data types must match the column headers of the input trace.
* `<redundancy>`: Specifies what kind of redundancy, if any, the unit has.  Must include a `type` attribute which can be either `"serial"` or `"parallel"` and a `count` attribute which must be a positive integer.
* `<trace>`: Indicates where to find the unit's trace file.  This is indicated by the `file` attribute and is a comma-separated file whose headers describe the operating point of the unit at each time step.  If a time step is missing a value, its value will be pulled from the corresponding `<default>` element. It also must include a comma-separated list of `failed` units. When a unit fails, each surviving unit switches its workload to the one described by the trace whose `failed` attribute contains the `name`s of the failed units. An empty string describes a system where no unit has failed.

The second section describes the *failure dependency graph* of the system. The failure dependency graph is a description of how failures propagate through the system from the architectural units up to the root. It consists of nested `<group>` and `<unit>` elements that track failures in their children. When enough failures occur in a group's children, the group itself fails, and when the top-level group fails, the entire system fails. This section might look something like this:
```
<group name="system" failures="0">
    <group name="a" failures="1">
        <unit name="core" />
        <unit name="accel" />
    </group>
    <group name="b" failures="0">
        <unit name="core" />
    </group>
</group>
```
A `<group>` element specifies a group of units or subgroups whose failures depend on each other. Like the `<unit>` element in the first section, it has a unique `name`. It also has a `failures` attribute, which specifies how many of its children can fail before it reports failure and must be nonnegative (0 means the group cannot tolerate failure). A `<group>` element may have any number of other `<group>` elements and any number of `<unit>` elements, all of which are its children. Unlike `<unit>` elements in the previous section, `<unit>` elements within a group only have a `name` attribute, and this attribute must correspond to one of the `unit`s in the first section.

Example configuration files can be found in the `example` directory.

### Trace Files
Each `<unit>` element in the first section of the configuration file must specify a trace file consisting of one or more time steps that indicate its operational point at each step. The first column of this file must be the time at which the step occurs, in arbitrary units. Each column after that describes an input value, such as supply voltage or clock frequency, at each time step. If a value is missing for a time step, OldSpot will use the default value specified by the corresponding `<default>` tag from the first section of the configuration file. An exampe trace might look like this:
```
time,activity,vdd,temperature,frequency,power
1,0.75,1.1,350,1000,1
```
The types of data required for each header depends on which aging mechanisms are selected to simulate, but generally include the following:
* activity: Activity of the unit, from 0 (never used) to 1 (used every cycle)
* vdd: Supply voltage in volts
* temperature: Temperature in Kelvins
* frequency: Clock frequency in MHz
* power: Power consumption in Watts
* current: Current drawn by the unit in Amperes
* current_density: Cross-sectional current density of the wires in the unit in Amperes/m<sup>2</sup>

Example trace files can be found in the `example` directory.

### Optional Parameters
OldSpot has several optional parameters that specify the number of Monte Carlo iterations to run, where to find technology and model parameters, which aging mechanisms to use, and what to output.  Run `./oldspot --help` for details on these parameters.

## Contact
If you have any questions about OldSpot or want to report bugs, please email Alec Roelke at <ar4jc@virginia.edu> or file an issue at the top of the page.
