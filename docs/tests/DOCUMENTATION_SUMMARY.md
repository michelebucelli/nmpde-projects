# NMPDE Projects: Documentation Summary

## Project Overview
The NMPDE (Numerical Methods for PDEs) repository implements space-time adaptivity algorithms for solving partial differential equations using the deal.II finite element library. The project features a comprehensive testing infrastructure with automated CI/CD pipelines, model validation, and Azure cloud deployment.

## Key Technologies
- **Language**: C++17 with deal.II library
- **Build System**: CMake (≥3.12)
- **CI/CD**: GitHub Actions (verification & scalability workflows)
- **Cloud**: Microsoft Azure with Terraform infrastructure-as-code
- **Testing**: CTest framework + Python validation (pyvista for VTU/PVTu analysis)
- **Containerization**: Docker (quay.io/pjbaioni/amsc_mk:2025)

## Testing Infrastructure
The project employs a multi-tiered testing approach:

1. **Verification Testing**: Automated build validation and unit tests via `./scripts/run_verification.sh`
2. **Model Validation**: Physical conservation law verification using Python validator analyzing VTU/PVTu mesh outputs
3. **Scalability Testing**: Thread-based performance benchmarking generating CSV result files
4. **Performance Monitoring**: System metric collection (CPU, memory, disk) with JSON/text reporting

## CI/CD Pipelines

**Verification Workflow** (automatic on push/PR):
- Builds project in Release mode with deal.II support
- Runs all CTest tests in parallel
- Downloads build artifacts for model validation stage
- Validates solution physical properties against expected reference values

**Scalability Workflow** (manual dispatch):
- Accepts optional custom executable parameter
- Runs benchmark tests across 1, 2, 4, 8, and 16 threads
- Generates scalability_results.csv with timing measurements
- Falls back to `ctest -R scalability` if no executable provided

## Model Validation Framework
- Supports CSV numeric files and VTU/PVTu visualization outputs
- Extracts `solution` field from mesh data using PyVista
- Computes volume-weighted mean verification against reference (default: 1/216 for unit cube)
- Configurable tolerance (default: 10⁻⁶) for relative error validation
- Ensures conservation laws are maintained by numerical scheme

## Azure Infrastructure
Terraform provisioning creates:
- **Resource Group**: Logical container (`nmpde-projects-rg`)
- **Container Registry (ACR)**: Image storage and GitHub Actions integration
- **Storage Account**: Three blob containers for reports, artifacts, and test results
- **Log Analytics**: Monitoring and diagnostic data aggregation

Service principal authentication enables GitHub Actions to deploy and manage Azure resources. All infrastructure is version-controlled and reproducible via Terraform configuration files.

## Quick Start Commands
```bash
# Local verification
./scripts/run_verification.sh

# Local model validation
./scripts/run_model_validation.sh build

# Local scalability testing
TEST_EXEC=build/space-time-adaptivity ./scripts/run_scalability.sh build

# Deploy Azure infrastructure
cd terraform && terraform apply

# Manual GitHub Actions workflow dispatch
# Open repository → Actions → Scalability Tests → Run workflow
```

## Key Documents
- **Full Documentation**: `Project_Documentation.tex` (5 pages) - comprehensive technical reference
- **CI/CD Guide**: `.github/workflows/` - workflow configurations
- **Azure Setup**: `docs/AZURE_SETUP.md` - detailed cloud deployment instructions
- **Validation Scripts**: `scripts/validate_model.py` - model correctness verification

---
*This documentation covers testing infrastructure, CI/CD automation, and cloud deployment for the NMPDE space-time adaptivity project.*
