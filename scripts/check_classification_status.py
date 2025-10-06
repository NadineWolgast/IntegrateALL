#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to check if classification was successful and determine if FusionCatcher needs to run.
Returns exit code 0 if FusionCatcher should run, 1 if classification was successful.
"""

import sys
import pandas as pd
import os

def needs_fusioncatcher(classification_file):
    """
    Check if FusionCatcher needs to run based on classification results.
    
    Args:
        classification_file: Path to final classification CSV
        
    Returns:
        bool: True if FusionCatcher should run (classification failed)
              False if classification was successful
    """
    
    if not os.path.exists(classification_file):
        print("Classification file not found: {}".format(classification_file))
        return True  # Run FusionCatcher if no classification exists
    
    try:
        # Read classification file
        df = pd.read_csv(classification_file)
        
        if df.empty:
            print("Classification file is empty")
            return True
            
        # Check WHO-HAEM5 column for successful classification
        if 'WHO-HAEM5' not in df.columns:
            print("WHO-HAEM5 column not found")
            return True
            
        who_haem5_value = df['WHO-HAEM5'].iloc[0] if len(df) > 0 else ""
        
        # Define exact patterns that indicate failed classification
        failed_patterns = [
            "IntegrateALL couldn't confirm the subtype",
            "IntegrateALL couldn't confirm the subtype in concordance with WHO or ICC classification"
        ]
        
        # Check for exact failed pattern matches or empty value
        who_haem5_str = str(who_haem5_value).strip()
        if who_haem5_str == "" or who_haem5_str in failed_patterns:
            print("Classification failed (WHO-HAEM5: {})".format(who_haem5_value))
            print("FusionCatcher will run to improve classification")
            return True
                
        # If we get here, classification was successful
        print("Classification successful (WHO-HAEM5: {})".format(who_haem5_value))
        print("FusionCatcher will be skipped")
        return False
        
    except Exception as e:
        print("Error reading classification file: {}".format(e))
        return True  # Run FusionCatcher on error


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_classification_status.py <classification_file>")
        sys.exit(2)
        
    classification_file = sys.argv[1]
    
    if needs_fusioncatcher(classification_file):
        print("RESULT: RUN_FUSIONCATCHER")
        sys.exit(0)  # Exit code 0 = run FusionCatcher
    else:
        print("RESULT: SKIP_FUSIONCATCHER") 
        sys.exit(1)  # Exit code 1 = skip FusionCatcher