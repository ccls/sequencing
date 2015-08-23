#!/bin/bash
su -l ec2-user -c 'nohup aws_1000genomes.sh --shutdown &> ~/aws_1000genomes.sh.log &'
