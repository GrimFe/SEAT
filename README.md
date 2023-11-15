# SEAT - Serpent EAsing Tool
<br>

![image](https://user-images.githubusercontent.com/69858810/207658811-29d7b404-f7ac-4d13-9979-f3ba2820afd2.png)

SEAT is a python package primarily meant to ease the user interface to the [Serpent](https://serpent.vtt.fi/mediawiki/index.php/Main_Page) Monte Carlo code inputs.
SEAT also comes with modules for processing Serpent results in the framwork of deterministic and stochastic uncertainty quantification.

## Structure
SEAT composes of several modules that are currently under development.
The advancement state of the modules ranges a lot.
Unit tests of SEAT are designed with [pytest](https://docs.pytest.org/en/7.4.x/) and we are working towards a for a full coverage.

SEAT composes of the following modules:
- `Serpent2InputParser` or [`sip`](https://github.com/GrimFe/SEAT/tree/master/SEAT/Serpent2InputParser): this module is not implemented yet; it would be meant to parse Serpent input files to `Serpent2InputWriter` objects;
- `Serpent2InputWriter` or [`siw`](https://github.com/GrimFe/SEAT/tree/master/SEAT/Serpent2InputWriter): the most developed SEAT module, designed to be a user friendly python interface to Serpent2 Input Files;
- `Serpent2UncertaintyPropagation` or [`sup`](https://github.com/GrimFe/SEAT/tree/master/SEAT/Serpent2UncertaintyPropagation): suite for the deterministic and stochastic uncertainty propagation with Serpent;

  ***

## 🔧 Installation
<br>
SEAT is not ready for pythoninc installation yet.
Temporarily, to use SEAT please clone [this](https://github.com/GrimFe/SEAT/tree/master) branch and import the modules you are using with:
```
import sys
sys.path.append(module_path)
```

## 📔 Documentation
<br>
We are working on a full documentation for SEAT, in the meanwhile, just refer to the object dcumentation in the form of docstrings.

## 🎮 Examples
<br>

## 📭 Contacts
<br>
* [**Federico Grimaldi**](https://github.com/GrimFe) - federico.grimaldi98@gmail.com
