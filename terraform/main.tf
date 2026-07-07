terraform {
  required_version = ">= 1.0"
  required_providers {
    azurerm = {
      source  = "hashicorp/azurerm"
      version = "~> 3.0"
    }
  }
}

provider "azurerm" {
  features {}
  subscription_id = var.subscription_id
}

# Resource Group
resource "azurerm_resource_group" "main" {
  name     = var.resource_group_name
  location = var.location
  tags = {
    Environment = "Development"
    Project     = "nmpde-projects"
    ManagedBy   = "Terraform"
  }
}

# Container Registry
resource "azurerm_container_registry" "acr" {
  name                = var.acr_name
  resource_group_name = azurerm_resource_group.main.name
  location            = azurerm_resource_group.main.location
  sku                 = "Basic"
  admin_enabled       = true
  tags = {
    Environment = "Development"
    Project     = "nmpde-projects"
  }
}

# Storage Account for performance reports and artifacts
resource "azurerm_storage_account" "storage" {
  name                     = var.storage_account_name
  resource_group_name      = azurerm_resource_group.main.name
  location                 = azurerm_resource_group.main.location
  account_tier             = "Standard"
  account_replication_type = "LRS"
  
  tags = {
    Environment = "Development"
    Project     = "nmpde-projects"
  }
}

# Blob container for performance reports
resource "azurerm_storage_container" "reports" {
  name                  = "performance-reports"
  storage_account_name  = azurerm_storage_account.storage.name
  container_access_type = "private"
}

# Blob container for build artifacts
resource "azurerm_storage_container" "artifacts" {
  name                  = "build-artifacts"
  storage_account_name  = azurerm_storage_account.storage.name
  container_access_type = "private"
}

# Blob container for test results
resource "azurerm_storage_container" "results" {
  name                  = "test-results"
  storage_account_name  = azurerm_storage_account.storage.name
  container_access_type = "private"
}

# Log Analytics Workspace for monitoring
resource "azurerm_log_analytics_workspace" "logs" {
  name                = "${var.resource_group_name}-logs"
  location            = azurerm_resource_group.main.location
  resource_group_name = azurerm_resource_group.main.name
  sku                 = "PerGB2018"
  retention_in_days   = 30
}

# Container App Environment for running MPI jobs
resource "azurerm_container_app_environment" "env" {
  name                       = "${var.resource_group_name}-env"
  location                   = azurerm_resource_group.main.location
  resource_group_name        = azurerm_resource_group.main.name
  log_analytics_workspace_id = azurerm_log_analytics_workspace.logs.id
}
