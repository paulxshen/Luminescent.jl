using VideoIO
a = zeros(UInt8, 350, 100)
open_video_out("a.mp4", a; framerate=24, encoder_options=(;)) do writer
    for t = 1:100
        write(writer, a)
    end
end
movie = VideoIO.openvideo("a.mp4", target_format=VideoIO.AV_PIX_FMT_GRAY8)
A = reinterpret(UInt8, read(movie))
@show size(A)
