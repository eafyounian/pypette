#!/bin/env julia

help = """
Tools for copy number analysis and visualization.

Usage:
  coverage tiled <bed_file> <grid_file>
  coverage cds <bam_file> <gtf_file>
  coverage logratio <test_wig> <ref_wig> <min_ref>

Options:
  -h --help       Show this screen.
  -s --step=N     Step size for window placement [default: window size / 2].
"""

using DocOpt, GZip

function ref_sequence_sizes(sam_path)
	ref_sizes = Dict{ASCIIString, Int64}()
	for line in eachline(`samtools view -H $(sam_path)`)
		m = match(r"@SQ\tSN:(\w+)\tLN:(\d+)", line)
		if m == nothing; continue; end
		ref_sizes[m.captures[1]] = int(m.captures[2])
	end
	return ref_sizes
end

function zopen(path)
	if path == "-"; return STDIN; end
	if ismatch(r"\.gz$", path)
		return GZip.open(path)
	else
		return open(path)
	end
end

info(msg) = println(STDERR, msg)
int(s::ASCIIString) = parse(Int, s)
int(s::SubString{ASCIIString}) = parse(Int, s)

type Chromosome
	name::ASCIIString
	first::Int32
	last::Int32
end

function coverage_tiled(bed_path, grid_path)
	# Read the karyotype, window size and step size from the grid file
	grid_file = zopen(grid_path)
	c = split(readline(grid_file), '\t')
	chr = Chromosome(c[1], int(c[2]) + 1, 0)
	chromosomes = [chr]
	winsize = int(c[3]) - int(c[2]); prev_pos = int(c[2])
	c = split(readline(grid_file), '\t')
	step = int(c[2]) - prev_pos; prev_pos = int(c[2])
	while !eof(grid_file)
		c = split(readline(grid_file), '\t')
		window = [int(c[2]), int(c[3])]
		@assert(window[2] - window[1] == winsize)
		if c[1] != chr.name
			chr.last = prev_pos
			chr = Chromosome(c[1], window[1] + 1, 0)
			push!(chromosomes, chr)
		else
			@assert(window[1] - prev_pos == step)
		end
		prev_pos = window[1]
	end
	chr.last = prev_pos
	win_overlap = winsize - step

	chr_list = join([replace(c.name, "chr", "") for c in chromosomes], ", ")
	info("Window size: $winsize")
	info("Step size: $step")

	bed_file = zopen(bed_path)
	cols = split(readline(bed_file), '\t')
	chr = filter(c -> c.name == cols[1], chromosomes)[1]
	cov = zeros(Int32, cld(chr.last - chr.first + 1, step))
	info("Analyzing chromosome $(chr.name)")
	while !eof(bed_file)
		cols = split(readline(bed_file), '\t')
		if cols[1] != chr.name
			for c in cov; println("$c"); end
			chr = filter(c -> c.name == cols[1], chromosomes)[1]
			cov = zeros(Int32, cld(chr.last - chr.first + 1, step))
			info("Analyzing chromosome $(chr.name)")
		end

		# Convert BED coordinates to 1-based inclusive genomic coordinates
		start = int(cols[2]) + 1
		stop = int(cols[3])

		# Convert coordinates to window indexes
		start = fld(start - chr.first - win_overlap, step) + 1
		start = max(start, 1)
		stop = cld(stop - chr.first, step) + 1
		#stop = min(stop, length(cov))   # FIXME: Unnecessary?

		cov[start:stop] += 1
	end
end

function coverage_cds(sam_path, gtf_path)
	chr_sizes = ref_sequence_sizes(sam_path)

	info("Constructing a map of coding regions...")
	coding = [chr => falses(size) for (chr, size) in chr_sizes]
	for line in eachline(zopen(gtf_path))
		if startswith(line, '#'); continue; end
		cols = split(line, '\t')
		if cols[3] != "CDS"; continue; end
		if length(cols[1]) > 5; continue; end
		ccds = get(coding, cols[1], nothing)
		if ccds == nothing; continue; end
		ccds[int(cols[4]):int(cols[5])] = true
	end

	info("Calculating a coverage histogram...")
	coverage = 0:200
	coverage_hist = zeros(length(coverage))
	chr = ""; pos = 0
	for line in eachline(`bedtools genomecov -d -split -ibam $(sam_path)`)
		cols = split(line, '\t')
		if cols[1] != chr
			chr = cols[1]
			pos = int(cols[2])
			ccds = coding[chr]
		else
			pos += 1
		end
		if ccds[pos]
			count = min(int(cols[3]), 200) + 1
			coverage_hist[count] += 1
		end
	end

	println("Coverage histogram:")
	println("===================")
	for (cov, count) in zip(coverage, coverage_hist)
		println("$(cov): $(count)")
	end
end


function coverage_logratio()
	for line in open(args["logratio"]["file"])

	end
end


args = docopt(help)
if args["tiled"]
	coverage_tiled(args["<bed_file>"], args["<grid_file>"])
elseif args["cds"]
	coverage_cds(args["<bam_file>"], args["<gtf_file>"])
elseif args["logratio"]
	coverage_logratio()
end
