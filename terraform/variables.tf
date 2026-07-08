variable "subscription_id" {
  description = "Azure Subscription ID"
  type        = string
  default     = "5d7fbbd4-db18-4472-915d-de671cc0f22a"
}

variable "resource_group_name" {
  description = "Name of the resource group"
  type        = string
  default     = "nmpde-projects-rg"
}

variable "location" {
  description = "Azure region for resources"
  type        = string
  default     = "italynorth"
}

variable "acr_name" {
  description = "Name of the Azure Container Registry (must be globally unique, alphanumeric only)"
  type        = string
  default     = "nmpdeprojectsacr"
}

variable "storage_account_name" {
  description = "Name of the storage account (must be globally unique, lowercase alphanumeric only, 3-24 chars)"
  type        = string
  default     = "nmpdeprojectsstorage"
}
