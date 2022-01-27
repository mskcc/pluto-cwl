#!/bin/bash
set -eux
dot fillout_diagram.dot -Tpdf -o fillout_diagram.pdf
dot samples_fillout_workflow.dot -Tpdf -o samples_fillout_workflow.pdf
