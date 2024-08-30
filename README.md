# Harmonic Flow Calculation and Power Quality Analysis Software

This repository contains the development of a software tool designed for analyzing and calculating harmonic flows and potential distortions in electrical systems. The project is divided into two distinct phases, each with specific objectives and advancements, but both focusing on power quality and harmonic distortion management.

## Phase 1: Development of a Program for Harmonic Flow Calculation

### Description

The first phase of the project focused on creating a software tool for calculating power flow and harmonic flow in electrical networks. The motivation for this phase was the increase in harmonic distortions due to the growing use of nonlinear devices in modern electrical grids.

**Objectives:**
- Develop a Python software for power flow calculation using the Newton-Raphson method.
- Implement harmonic flow calculation to provide a comprehensive analysis of operation in both fundamental and harmonic frequencies.
- Validate the results obtained by comparing them with data from an established software in the electrical sector.
- Lay the groundwork for studying and analyzing harmonic responsibility sharing among different agents and consumers.

**Keywords:** Harmonic sharing, Software development, Power flow, Harmonic flow, Power quality.

### Usage

1. **Prerequisites:** Python 3.12.5 (latest) and math related libraries.
2. **Installation:** Clone this repository and install the required dependencies listed in the `requirements.txt` file.

    ```bash
    pip install -r requirements.txt
    ```

3. **Execution:** Use the scripts in the `main.py` file to calculate power and harmonic flow. Refer to the README within that directory for specific usage details.

## Phase 2: Implementation of Blind Source Separation Algorithm (FAST ICA)

### Description

The second phase of the project extended the functionality of the previously developed software by incorporating the Blind Source Separation (BSS) technique using the FAST ICA algorithm. This enhancement allows for a more detailed analysis of each source’s individual contribution to harmonic distortion.

**Objectives:**
- Implement the FAST ICA algorithm for blind source separation and analysis of individual harmonic contributions.
- Enhance the software’s capabilities to provide an innovative and non-invasive tool for managing harmonic distortions.
- Conduct bibliographic studies and specific case analyses to validate the effectiveness of the new feature.

**Keywords:** Blind source separation, FAST ICA, Harmonic distortion, Power quality, Software development.

### Usage

1. **Prerequisites:** Python 3.12.5 (latest),  math related, signal processing and analysis libraries (see `requirements.txt`).
2. **Installation:** Clone this repository and install the required dependencies listed in the `requirements.txt` file.

    ```bash
    pip install -r requirements.txt
    ```

3. **Execution:** Use the scripts in the `main.py` file to apply the FAST ICA algorithm and analyze harmonic contributions. Refer to the README within that directory for specific usage details.

## Contributing

If you wish to contribute to the project, please follow these steps:

1. Fork this repository.
2. Create a new branch for your changes (`git checkout -b my-contribution`).
3. Commit your changes (`git commit -am 'Add new feature'`).
4. Push your branch to the repository (`git push origin my-contribution`).
5. Open a pull request with a clear description of your changes.

## License

This project is licensed under the [MIT License](LICENSE). See the LICENSE file for more details.

## Contact

For more information about the project or for inquiries, please contact:

- **Breno Marques Freitas** - [brenomarques137@gmail.com](mailto:your-email@example.com)
- **Andréia Crico dos Santos** - [andreiacrico@iftm.edu.br](mailto:your-email@example.com)

Thank you for visiting the repository!
