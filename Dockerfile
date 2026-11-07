FROM docker.io/mambaorg/micromamba:1.5.10

# SciLifeLab Serve requires a non-root user with UID 1000
USER root
RUN usermod -u 1000 $MAMBA_USER && groupmod -g 1000 $MAMBA_USER \
    && chown -R 1000:1000 /opt/conda /home/$MAMBA_USER
USER $MAMBA_USER

# Install Python, RDKit, Shiny etc. from environment.yml
COPY --chown=1000:1000 environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# The app and its startup script
WORKDIR /app
COPY --chown=1000:1000 chase/ chase/
COPY --chown=1000:1000 start-script.sh start-script.sh
RUN chmod +x start-script.sh

EXPOSE 8765
CMD ["./start-script.sh"]
