#!/bin/bash
current_time=$(date | tr -s " " "\t" | cut -f 4 | cut -d ":" -f 1,2)
echo "The time is $current_time.
I'm glad to see you're making good use of it :)"
