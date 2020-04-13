#!/usr/bin/env awk -f

/^[^#]/ { 
    acov = ($4 - $5 - $6) / $4 * 100;
    match($13, "^(.-)+");
    gapl = RLENGTH;
    match($13, "(.-)+$");
    gapl = RLENGTH;
    print $0, acov, gapl, gapr;
}
