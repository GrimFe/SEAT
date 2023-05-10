from dataclasses import dataclass

__author__ = "Federico Grimaldi"
__all__ = [
    "Detector",
    ]

@dataclass(slots=True)
class Detector:
    name: str
