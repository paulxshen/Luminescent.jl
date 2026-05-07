using VideoIO
T = UInt8
target_pix_fmt = VideoIO.AV_PIX_FMT_GRAY8
MAX = typemax(T)
@time open_video_out("_.mp4", zeros(T, 100, 100); target_pix_fmt) do writer
    for i = vcat(1:50, 1:50)
        write(writer, convert.(T, round.(i / 100 * MAX * ones(100, 100))))
    end
end

a = []
movie = VideoIO.openvideo("_.mp4")
vr = VideoIO.seekstart(movie)
for i = 1:100
    global x = read(vr)
    y = reinterpret(T, x) / MAX
    @show maximum(y) - minimum(y)
    @show y[1]
    push!(a, y)
end