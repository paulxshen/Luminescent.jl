using VideoIO

# Construct a AVInput object to access the video and audio streams in a video container
# io = VideoIO.open(video_file)
io = VideoIO.testvideo("annie_oakley") # for testing purposes

# Access the video stream in an AVInput, and return a VideoReader object:
f = VideoIO.openvideo(io) # you can also use a file name, instead of a AVInput
f = reverse(f)
img = read(f)

while !eof(f)
    read!(f, img)
    # Do something with frames
end
close(f)