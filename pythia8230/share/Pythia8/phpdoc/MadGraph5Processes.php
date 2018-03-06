<html>
<head>
<title>MadGraph5 Processes</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='MadGraph5Processes.php'>
 
<h2>MadGraph5 Processes</h2> 
 
Here we will describe two special ways to make use of MadGraph5 and 
MadGraph5_aMC@NLO [<a href="Bibliography.php#refAlw11" target="page">Alw11</a>,<a href="Bibliography.php#refAlw14" target="page">Alw14</a>]inside PYTHIA, either by exporting 
Madgraph process code or by wrapping the MadGraph5_aMC@NLO generator as 
a PYTHIA Les Houches interface. 
 
<p/> 
Of course, MadGraph5 can also output files of parton-level events 
according to the <?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>LHEF</a> standard, 
that can be read in and processed further by PYTHIA 8. This is the 
most commonly used approach, and requires no further description here. 
 
<h3>MadGraph5 code inside PYTHIA</h3> 
 
By far the easiest way to implement new processes into PYTHIA 8 is 
by using the matrix-element generator MadGraph5. This program has 
an option to output the results of a matrix-element calculation 
as a set of PYTHIA 8 C++ classes (plus further auxiliary code), 
that can then be linked and used as 
<?php $filepath = $_GET["filepath"];
echo "<a href='SemiInternalProcesses.php?filepath=".$filepath."' target='page'>";?>semi-internal</a> processes, 
meaning they are handled identically with normal internal ones. 
This way, MadGraph5 can be used to implement processes from 
any model that can be written in  terms of a Lagrangian. Any 
<i>2 &rarr; 1</i>, <i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes 
can be implemented, the limit being set by the absence of efficient 
phase space generator algorithms for higher multiplicities in 
PYTHIA. Features such as <i>s</i>-channel resonances are 
automatically implemented in the process classes.  Besides the process 
library and necessary model files, also an example main program is 
generated for each set of processes, which can be easily modified to 
perform the desired analyses. 
 
<p/> 
In order to create a PYTHIA 8 process library with MadGraph5, first 
download the MadGraph5 package from 
<a href="https://launchpad.net/madgraph5" target="page"> 
https://launchpad.net/madgraph5</a>, and untar the package. You can 
then specify the location of your <code>pythia81xx</code> directory 
in the file <code>input/mg5_configuration.txt</code>: 
<br/><code>pythia8_path = ./pythia81xx</code> 
<br/>The location can be either relative (to the directory 
<code>MadGraph5_v_x_x_x/.</code>) or absolute. 
 
<p/> 
For any model that is already implemented in the MadGraph5 package, 
you can directly use the model. Start the MadGraph5 interface 
<code>bin/mg5</code>, and do: 
<pre> 
import model model_name 
generate your_process_in_mg5_syntax 
add process your_next_process_in_mg5_syntax 
... 
output pythia8 [path_to_pythia81xx_directory] 
</pre> 
 
<p/> 
For examples of MG5 process syntax, please see 
<a href="http://madgraph.phys.ucl.ac.be/EXAMPLES/example_mg5.html" 
target="page">http://madgraph.phys.ucl.ac.be/EXAMPLES/example_mg5.html</a> 
or type <code>help generate</code>. If you specified the path to the 
<code>pythia81xx</code> directory in the <code>mg5_configuration</code> 
file, you do not need to enter it in the <code>output</code> command. 
 
<p/> 
If your preferred model is found on the FeynRules model wiki page, 
<a href="http://feynrules.irmp.ucl.ac.be/wiki/ModelDatabaseMainPage" 
target="page">http://feynrules.irmp.ucl.ac.be/wiki/ModelDatabaseMainPage</a>, 
download the UFO (Universal FeynRules Output) tar file for the model, 
untar in the <code>models/</code> directory, and use as above. 
 
<p/> 
If you want to implement a new model which has not yet been implemented, 
you can do this either using the Mathematica package FeynRules (see 
<a href="http://feynrules.irmp.ucl.ac.be/" target="page"> 
http://feynrules.irmp.ucl.ac.be/</a>) or directly edit the UFO model 
files of the most similar model in the <code>models/</code> directory. 
 
<p/> 
The resulting output from the <code>output pythia8</code> command is: 
<ul> 
<li>A process directory <code>Processes_modelname</code> with the 
model information and the files needed for all processes defined for 
this model, placed in the <code>pythia81xx</code> main directory. 
The model files are <code>Parameters_modelname.h/cc</code> and 
<code>HelAmps_modelname.h/cc</code>, and the process files for each 
process class (with the same mass, spin and color of the initial/final 
state particles) are called <code>Sigma_modelname_processname.h/cc</code>. 
The directory also contains a <code>makefile</code> and a model parameter 
file <code>param_card_modelname.dat</code>.</li> 
<li>An example main program in the directory <code>examples/</code> 
(in the <code>pythia81xx</code> main directory) called 
<code>main_modelname_N.cc</code> and a corresponding makefile 
<code>Makefile_modelname_N</code>. This main program links in the 
process classes in the process directory described above. To run the 
example main program, just go to the <code>examples/</code> 
directory and run 
<br/><code>make -f Makefile_modelname_N</code> 
<br/>or run <code>launch</code> directly inside the MadGraph5 
command line interface.</li> 
</ul> 
 
<p/> 
Note that in order for PYTHIA 8 to be able to automatically decay any 
new particles, it is necessary to specify the branching ratios of the 
particles in the <code>param_card</code> file, see 
[<a href="Bibliography.php#refSka04" target="page">Ska04</a>,<a href="Bibliography.php#refAlw07" target="page">Alw07</a>] for details. 
 
<p/> 
For further technical details, please see the MadGraph5 release paper 
[<a href="Bibliography.php#refAlw11" target="page">Alw11</a>] and the 
<?php $filepath = $_GET["filepath"];
echo "<a href='SemiInternalProcesses.php?filepath=".$filepath."' target='page'>";?>semi-internal</a> processes page. 
 
<p/> 
Currently the standard way of interfacing is to use the LHEF standard 
with an intermediate event file. The advantage is that then the 
MadGraph5 phase space generator can be used, which opens up for 
processes with more than three particles in the final state. The 
disadvantages are that it is less easy to mix and match with existing 
PYTHIA processes, and that one needs to regenerate and store large LHEF 
files for different  kinematics cuts or parameter values. 
 
<p/> 
Please cite the MadGraph5 release paper [<a href="Bibliography.php#refAlw11" target="page">Alw11</a>] if you use 
MadGraph5 to generate process libraries for PYTHIA 8. 
 
<h3>MadGraph5_aMC@NLO executable inside PYTHIA</h3> 
 
The <code>Pythia::setLHAupPtr(LHAup* lhaUpPtr)</code> method allows 
a Pythia generator to accept a pointer to an object derived from the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>LHAup</a></code> base class. 
Such an object will be initialized from within Pythia, and be called 
repeatedly to generate the next parton-level event, using the LHA 
specification as a standard to transfer the relevant information back 
to Pythia. Properly constructed, the operation of an <code>LHAup</code> 
object thus is almost completely hidden from the user, and generates 
events almost like an ordinary internal Pythia process. 
 
<p/> 
The <code>LHAupMadgraph</code> is precisely such a class, derived from 
<code>LHAup</code>, that contains the code needed to wrap a 
MadGraph5_aMC@NLO executable. Thereby the generation of Madgraph 
processes from within Pythia becomes straightforward. An explicit 
example is provided in <code>main34.cc</code>. We describe some of the 
key elements used there and in the general case. 
 
<a name="anchor1"></a>
<p/><strong>LHAupMadgraph::LHAupMadgraph(Pythia* pythia, bool match = true, string dir = &quot;madgraphrun&quot;, string exe = &quot;mg5_aMC&quot;) &nbsp;</strong> <br/>
creates an instance of the <code>LHAupMadgraph</code> class. 
<br/><code>argument</code><strong> pythia </strong>  :  pointer to the <code>Pythia</code> instance, 
such that some of its facilities can be used inside the interface. 
   
<br/><code>argument</code><strong> match </strong> (<code>default = <strong>on</strong></code>) :  should be true if jet matching is 
requested. For tree-level generation MLM matching is used, while 
FxFx matching is used for aMC@NLO generation. This is set up in 
<code>LHAupMadgraph::setInit()</code>, which could be modified to 
represent other matching strategies or parameter values. 
   
<br/><code>argument</code><strong> dir </strong> (<code>default = <strong>madgraphrun</strong></code>) :  the name of the run 
directory, into which MadGraph puts its (intermediate) results. 
   
<br/><code>argument</code><strong> exe </strong> (<code>default = <strong>mg5_aMC</strong></code>) :  the name of the MadGraph5_aMC@NLO 
executable that <code>LHAupMadgraph</code> is meant to wrap. In additon 
it may be necessary to prepend the full pathname of the executable: 
<code>"(something)/MG5_aMC_v2_3_3/bin/mg5_aMC"</code>. 
   
   
 
<a name="anchor2"></a>
<p/><strong>bool LHAupMadgraph::readString(string line, Stage stage = Auto) &nbsp;</strong> <br/>
allows the user to send commands to MadGraph. 
<br/><code>argument</code><strong> line </strong>  :  the command to be sent to MadGraph. 
Any string begining with <code>"configure "</code> is used for the initial 
MadGraph configuration with <code>"configure "</code> stripped from the 
begining. In general, only the process and run settings need to be 
provided. Run settings must begin with <code>" set"</code>; note the 
leading space. The output and launch commands, random seed, and shower 
choice are automatically handled. For example, the following will produce 
di-muon events from 13 TeV proton proton collisions at NLO in QCD: 
<br/><code>readString("generate p p > mu+ mu- [QCD]");</code> 
   
<br/><code>argument</code><strong> stage </strong> (<code>default = <strong>Auto</strong></code>) :  if the stage is set to Auto, 
commands beginning with <code>" set"</code> are used in the launch 
stage, commands begining with <code>"configure"</code> are used in the 
configuration stage, and all remaining commands 
(excluding <code>output</code> and 
<code>launch</code>) are used in the generate stage. Output, launch, seed, 
and shower commands are automatically handled. If the user wishes to 
override commands, then the stage can be specified. This will prevent any 
automatically generated commands from being used for that stage. This 
should only be done if the user understands what additional commands are 
needed. 
   
   
 
<a name="anchor3"></a>
<p/><strong>void LHAupMadgraph::setEvents(int events) &nbsp;</strong> <br/>
the number of events to generate per MadGraph run. Normally does not 
need to be set, but defaults to 10000. 
   
 
<a name="anchor4"></a>
<p/><strong>void LHAupMadgraph::setSeed(int seed, int runs = 30081) &nbsp;</strong> <br/>
the random seed (sequence), normally not needed to be set explicitly. 
If the random seed is negative (default of -1), then the MadGraph 
seed is taken as the Pythia parameter <code>"Random:seed"</code>, which 
must be greater than 0. If the maximum number of allowed runs is exceeded 
(default of 30081) an error is thrown. The seed for a MadGraph run is set as: 
<br/> (random seed - 1) * (maximum runs) + (number of runs) + 1. 
<br/>MadGraph can only handle random seeds up to 30081 * 30081. So, with 
this strategy, one can generate Pythia jobs with seeds from 1 to 30081, 
with each job running MadGraph less than 30081 times, and ensure a fully 
statistically independent sample. If more than 30081 jobs are needed, then 
the maximum allowed runs can be lowered accordingly, and if need be, 
setEvents can be used to increase the number of events generated per run. 
   
 
<a name="anchor5"></a>
<p/><strong>void LHAupMadgraph::setJets(int jets) &nbsp;</strong> <br/>
Set the number maximum number of jets generated by MadGraph. 
If negative (default of -1) then the number of jets is determined 
automatically, to be the maximum number of jets produced at leading order. 
   
 
<p/> 
Note that GZIP support must be enabled in the Pythia executable, so use 
the <code>--with-gzip</code> option in the <code>configure</code> step 
before Pythia compilation. 
 
<p/> 
Events are generated with MadGraph utilizing the 
<a href="https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/GridDevelopment" 
target="_top">gridpack</a> method for MadGraph5 and an 
<a href="https://answers.launchpad.net/mg5amcnlo/+question/243268" 
target="_top">equivalent method</a> for aMC@NLO. Consequently the run 
directory, <code>"madgraphrun"</code> by default, does not need to be deleted 
between independent runs with the same configuration (excluding random 
seeds). Indeed, keeping the directory significantly speeds the generation 
process, particularly for NLO generation with aMC@NLO as the grid 
initialization can be skipped after the initial run. 
 
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 