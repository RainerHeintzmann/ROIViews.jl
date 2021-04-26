using IndexFunArrays, ROIViews

# a random position within the size limit sz
randpos(sz, dims) = Tuple((d <= length(sz)) ? rand(1:sz[d]) : 1 for d in 1:dims)

@testset "ROIView" begin
    for d in 1:5
        sz = Tuple(9 .+ ones(Int,d));
        data = idx(sz,offset=CtrCorner);
        N = 3
        for s = 1:d # smaller ROI dimensions. s is ROI dimensionality
            ROIpos = Tuple(Tuple(rand(1:10, s)) for n in 1:N) 
            ROIsize = Tuple(rand(1:10,s))
            pad_val = 0 .* data[1];
            rois = ROIView(data, ROIpos, ROIsize, pad_val=pad_val);
            for p = 1:N # select a ROI number
                p_roi = randpos(ROIsize,d) # position in the ROI data
                @test size(rois) == (expand_size(ROIsize,size(data))...,N)
                roi_val = rois[p_roi...,p]
                pos = p_roi .+ Tuple((n<=s) ? ROIpos[p][n] .- ROIsize[n].÷2 : 0 for n=1:d)
                if Base.checkbounds(Bool, data, pos...)
                    data_val = data[pos...]
                    @test roi_val .+ 1 == pos
                else
                    data_val = pad_val
                end
                @test data_val == roi_val
            end
        end
    end
end
