#!/usr/bin/env awk -f
/^>/{ name=$0 } 
!/(^>)|N/{
	print name
	print $0
}
