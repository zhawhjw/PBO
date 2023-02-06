# Preferential Bayesian Optimizer

Customized C++ PBO version based on [Sequential-Line-Search](https://github.com/yuki-koyama/sequential-line-search). Yuki Koyama (<yuki@koyama.xyz>) is the original author.

### Prerequisites

#### MacOS
```
brew install cmake eigen
```

#### Linux
```
sudo apt install cmake libeigen3-dev
```

#### Windows 10
Install visual studio 2019 (other versions are not tested) with the MSVC c++ compiler.
The 'eigen-3.4.0' folder is used for Windows 10 only for build convenience.

## How to Compile and Run (C++)

#### MacOS or Linux
Use [cmake](https://cmake.org/) (3.3 or above) and then
```bash
git clone https://github.com/zhawhjw/PBO.git
cd PBO
mkdir build
cd build
cmake ../
make
```
Stay at the 'build' directory, then you can run the applications by (there is a hard-coded relative path in main.c, so make sure you at the correct directory which is achieved by following commands),
```
cd demos/sequential_line_search_nd/
./SequentialLineSearchNd --csv ABSOLUTE_PATH_TO_CSV
```
#### Windows 10
Use Microsoft Visual Studio 2019 to open the project folder, the MS should automatically build the project and generate the "out" folder in project folder.
Navigate to the location of built object file
```
cd out/build/your-debugger (usually it is x64-Debug)/demos/sequential_line_search_nd
./SequentialLineSearchNd.exe --csv ABSOLUTE_PATH_TO_CSV
```

#### --csv ABSOLUTE_PATH_TO_CSV
A parameter passed to built object file. The default value is 
```
ABSOLUTE_PATH_TO_CSV = PATH/Corridor_Scenario.csv
```
Or you can choose to type any of below ones:

```
PATH/Corridor_Scenario.csv
PATH/Hallway_Scenario.csv
PATH/Quad_Scenario.csv
```

### Data
The data is placed at folder "csv". Each row represents one parameter combination.
For columns: the last one is the preference percentage and all previous columns are parameters.
The preference tells how many users favor the current experiment. Be notice that there are multiple rows to be put together to become an experiment.
The preference percentage should be added up to 100% in one experiment. Below is how we organized our experiments, the digit represents the index of the row in csv file:
```c++
vector<int> first = {0, 1} ;
vector<int> second = {2, 3, 4};
vector<int> third = {5, 6, 7};
vector<int> forth = {8, 9, 10};
vector<int> fifth = {11, 12, 13};
vector<int> sixth = {14, 15};
```

[//]: # (### Source Code)

[//]: # ()
[//]: # (The main source code at directory "demos/sequential_line_search_nd/main.c")

[//]: # ()
[//]: # (Line 154 of "main.c" contains the path to different scenario)

[//]: # ()
[//]: # (```)

[//]: # (std::string strPathCSVFile = "../../../csv/Corridor_Scenario.csv";)

[//]: # (```)

## License

MIT License.


