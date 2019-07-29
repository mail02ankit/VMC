#This script add one blank line after one complete iteration.
START{temp=x;}{if(temp!=$1){{print "";}{temp=$1;}};}{print ;}
