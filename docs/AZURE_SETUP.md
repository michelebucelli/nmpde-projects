# Azure CI/CD Setup Guide

This guide will help you set up GitHub Actions to run your C++ project on Azure with MPI support and performance reporting.

## 🏗️ Infrastructure Overview

The Terraform setup provisions:
- **Azure Container Registry (ACR)**: Mirrors your quay.io container
- **Storage Account**: Stores performance reports, artifacts, and test results
- **Container App Environment**: Runs MPI workloads
- **Log Analytics**: Monitors job execution

## 📋 Prerequisites

1. **Azure Account** with an active subscription
2. **Azure CLI** installed locally (`brew install azure-cli` on macOS)
3. **Terraform** installed (`brew install terraform`)
4. **GitHub repository** with admin access

## 🚀 Setup Steps

### Step 1: Azure Service Principal Setup

Create a service principal for GitHub Actions to authenticate with Azure:

```bash
# Login to Azure
az login

# Set your subscription
az account set --subscription 5d7fbbd4-db18-4472-915d-de671cc0f22a

# Create service principal
az ad sp create-for-rbac \
  --name "github-actions-nmpde" \
  --role contributor \
  --scopes /subscriptions/5d7fbbd4-db18-4472-915d-de671cc0f22a \
  --sdk-auth
```

This will output JSON credentials. **Save this output** - you'll need it for GitHub Secrets.

### Step 2: Configure GitHub Secrets

Go to your GitHub repository settings → Secrets and variables → Actions → New repository secret

Add the following secret:

**Name:** `AZURE_CREDENTIALS`  
**Value:** Paste the entire JSON output from Step 1

The JSON should look like:
```json
{
  "clientId": "...",
  "clientSecret": "...",
  "subscriptionId": "5d7fbbd4-db18-4472-915d-de671cc0f22a",
  "tenantId": "...",
  "activeDirectoryEndpointUrl": "...",
  "resourceManagerEndpointUrl": "...",
  "activeDirectoryGraphResourceId": "...",
  "sqlManagementEndpointUrl": "...",
  "galleryEndpointUrl": "...",
  "managementEndpointUrl": "..."
}
```

### Step 3: Deploy Infrastructure with Terraform

From your local machine:

```bash
# Navigate to terraform directory
cd terraform

# Initialize Terraform
terraform init

# Preview changes
terraform plan

# Apply infrastructure
terraform apply
```

Enter `yes` when prompted. This will create:
- Resource Group: `nmpde-projects-rg`
- Container Registry: `nmpdeprojectsacr`
- Storage Account: `nmpdeprojectsstorage`

### Step 4: Verify Infrastructure

```bash
# Check resource group
az group show --name nmpde-projects-rg

# List all resources
az resource list --resource-group nmpde-projects-rg --output table

# Verify ACR
az acr list --resource-group nmpde-projects-rg --output table

# Verify Storage
az storage account list --resource-group nmpde-projects-rg --output table
```

### Step 5: Trigger Your First Azure Workflow

Push to your repository or manually trigger the workflow:

```bash
# Push changes
git add .
git commit -m "Add Azure CI/CD infrastructure"
git push

# Or trigger manually via GitHub UI:
# Go to Actions → Azure MPI Performance Tests → Run workflow
```

## 📊 Workflow Features

### Automated Tests
- **Build & Test**: Compiles project with CMake and runs unit tests
- **MPI Performance**: Runs with configurable process counts (2, 4, 8)
- **Performance Metrics**: Collects CPU, memory, disk, timing data
- **Artifact Storage**: Uploads reports to Azure Storage

### Manual Triggers
You can trigger the workflow manually with custom parameters:
- Number of MPI processes
- Test type (convergence, scalability, full)

### Scheduled Runs
The workflow runs automatically every night at 2 AM UTC.

## 📂 Accessing Results

### Via GitHub Actions UI
1. Go to Actions tab in GitHub
2. Select a workflow run
3. View job summaries and download artifacts

### Via Azure Storage
```bash
# List all performance reports
az storage blob list \
  --account-name nmpdeprojectsstorage \
  --container-name performance-reports \
  --auth-mode login \
  --output table

# Download a specific report
az storage blob download \
  --account-name nmpdeprojectsstorage \
  --container-name performance-reports \
  --name run-123/performance_20260222_120000.json \
  --file ./local_report.json \
  --auth-mode login
```

### Via Azure Portal
1. Go to [Azure Portal](https://portal.azure.com)
2. Navigate to `nmpde-projects-rg` resource group
3. Open `nmpdeprojectsstorage` storage account
4. Browse `performance-reports` container

## 🔧 Configuration

### Modify Infrastructure
Edit [`terraform/variables.tf`](terraform/variables.tf) to change:
- Resource names
- Regions
- SKUs/pricing tiers

Then run:
```bash
cd terraform
terraform plan
terraform apply
```

### Modify Workflow
Edit [`.github/workflows/azure-deploy.yml`](.github/workflows/azure-deploy.yml) to:
- Change MPI process counts
- Add custom build flags
- Modify test parameters
- Add new performance metrics

### Adjust MPI Processes
Edit the workflow file line with `mpi_processes` default value, or pass it as input when manually triggering.

## 🧪 Local Testing

Test the performance script locally:

```bash
# Build your project first
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Run performance tests
chmod +x scripts/run_performance_tests.sh
scripts/run_performance_tests.sh build reports
```

## 📈 Performance Reports

Reports include:
- **Timestamp**: When the test ran
- **System Info**: CPU, memory, disk specs
- **Duration**: Execution time per test
- **Memory Usage**: Peak memory consumption
- **Exit Codes**: Success/failure status
- **Test Output**: Full console output

### Example Report Structure
```json
{
  "timestamp": "2026-02-22T10:30:00Z",
  "hostname": "azure-runner",
  "cpu_count": 8,
  "memory_gb": 16,
  "tests": [
    {
      "name": "convergence",
      "duration_seconds": 45.2,
      "max_memory_kb": 512000,
      "exit_code": 0
    }
  ]
}
```

## 🐛 Troubleshooting

### Authentication Errors
```bash
# Verify service principal
az ad sp show --id <clientId-from-secret>

# Check role assignments
az role assignment list --assignee <clientId-from-secret>
```

### Container Registry Issues
```bash
# Test ACR access
az acr login --name nmpdeprojectsacr

# Verify imported image
az acr repository list --name nmpdeprojectsacr
```

### Storage Access Issues
```bash
# Check storage account
az storage account show --name nmpdeprojectsstorage

# Verify containers
az storage container list --account-name nmpdeprojectsstorage --auth-mode login
```

### Workflow Failures
1. Check the GitHub Actions logs for detailed error messages
2. Verify all secrets are set correctly
3. Ensure Azure resources exist
4. Check service principal permissions

## 🧹 Cleanup

To remove all Azure resources:

```bash
cd terraform
terraform destroy
```

Or delete the entire resource group:
```bash
az group delete --name nmpde-projects-rg --yes
```

## 🔐 Security Best Practices

- ✅ Never commit `AZURE_CREDENTIALS` to git
- ✅ Use managed identities when possible
- ✅ Regularly rotate service principal credentials
- ✅ Set minimum required permissions (principle of least privilege)
- ✅ Enable Azure Storage encryption (enabled by default)
- ✅ Use Azure Key Vault for sensitive data (optional upgrade)

## 📚 Additional Resources

- [Azure Container Registry Documentation](https://docs.microsoft.com/azure/container-registry/)
- [GitHub Actions for Azure](https://github.com/Azure/actions)
- [Terraform Azure Provider](https://registry.terraform.io/providers/hashicorp/azurerm/latest/docs)
- [Azure Storage Documentation](https://docs.microsoft.com/azure/storage/)

## 🆘 Support

For issues:
1. Check workflow logs in GitHub Actions
2. Review Terraform state: `terraform show`
3. Check Azure Portal for resource status
4. Review container logs in Azure

## 🎯 Next Steps

1. **CUDA Support**: Add GPU-enabled Container Instances for CUDA tests
2. **Custom Metrics**: Extend performance scripts with domain-specific metrics
3. **Dashboards**: Create Azure Dashboard for visualization
4. **Alerts**: Set up Azure Monitor alerts for failures
5. **Cost Optimization**: Review Azure Advisor recommendations
