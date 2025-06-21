# GenFitExamples

This repository provides example code and documentation for working with [GenFit2](https://github.com/GenFit/GenFit), a generic track fitting toolkit widely used in particle physics for detector reconstruction.
The examples are intended to help new users get started with GenFit2 and demonstrate common use cases, such as calculating fit residuals and integrating GenFit2 into analysis workflows.

## What is GenFit2?

[GenFit2](https://github.com/GenFit/GenFit) is a C++ library designed for generic track fitting.
It supports a variety of detector geometries and material effects, making it a flexible tool for high-energy physics experiments.
GenFit2 is open source and can be integrated into larger frameworks or used as a standalone toolkit for prototyping and educational purposes.

## Repository Structure

- `GFResidualExample/`: Contains an example demonstrating how to calculate and analyze fit residuals using GenFit2.
- Additional folders may provide further examples for other GenFit2 features as the repository evolves.

## Getting Started

1. **Clone the Repository**
   ```bash
   git clone https://github.com/pbielefeldt/GenFitExamples.git
   cd GenFitExamples
   ```

2. **Install GenFit2**  
   Follow the [GenFit2 installation instructions](https://github.com/GenFit/GenFit#installation) to build and install the library.
   Make sure to set up your environment variables so that GenFit2 headers and libraries are discoverable.

3. **Build the Examples**  
   Follow the instructions in each example's README.

## Usage Notes

- These examples assume you have a working installation of GenFit2 and required dependencies (such as ROOT).
- Source code is provided as a reference and starting point; adapt it as needed for your specific use case or detector setup.
- For further details and advanced usage, see the [GenFit2 documentation](https://github.com/GenFit/GenFit/wiki).


---
**Acknowledgments:**  
This project is not affiliated with the official GenFit2 developers.
