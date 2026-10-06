"""Content-addressed provenance for trained models and resource-based methods."""
from .registry import describe_model, write_provenance

__all__ = ["describe_model", "write_provenance"]
