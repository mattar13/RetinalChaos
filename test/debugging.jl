using PyPlot

# Generate the data for the animation
data = rand(64, 64, 12000)

# Set up the figure and axes for the animation
fig, ax = subplots()

# Initialize the image plot
im = imshow(data[:,:,1])

# Function to update the plot for each frame of the animation
function update(i)
     println(i)
     im.set_data(data[:,:,i+1])
     return im
end

# Set up the animation using PyPlot's animation module
animation = anim.FuncAnimation(fig, update, frames=12000, interval=50)

# Save the animation as an MP4 file using FFmpeg
animation.save("animation.mp4", dpi=80, writer="ffmpeg")